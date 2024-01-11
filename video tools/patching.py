# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 12:24:59 2022

@author: Yuchen Wang
"""

import argparse
import os
import itertools
from shutil import copyfile

parser = argparse.ArgumentParser()
parser.add_argument('config', type=str, help='the configure file')
parser.add_argument('-i', '--info', dest='showonly', action='store_true', help='do not perform any commands, only show them')
parser.add_argument('-e', '--temp-ext', dest='tempext', default='mp4', help='extension for temp files, e.g. "mp4", "flv" (no "."!)')
parser.add_argument('-c:a', dest='ca', default='copy')
parser.add_argument('-c:v', dest='cv', default='copy')
parser.add_argument('-crf', dest='crf', default=10, help='constant rate factor (crf) for h264.')
args = parser.parse_args()

inf = 99999

def parse_time(s):
    if type(s) is str and':' in s:
        ts = s.split(':')
        secs = sum([float(t) * 60**i for i, t in enumerate(ts[::-1])])
        return secs
    elif s in ['end', 'END']:
        return inf
    else:
        return float(s)
    
def run(cmd):
    print(f'[RUN] {cmd}')
    if not args.showonly: os.system(cmd)

def copy(src, dst):
    print(f'[COPY] {src} -> {dst}')
    if not args.showonly: copyfile(src, dst)

input_videos = {}
final_cmds = []

outname = None

with open(args.config, encoding='utf8') as f:
    for line in f.readlines():
        line = line.replace('\n', '')
        match line.split('\t'):
            case ('#IN', inname, *_):
                input_videos[inname] = []
            case ('*clip', *clipinfo):
                input_videos[inname].append(clipinfo)
            case (('*rm' | '*remove'), rmtime, *_):
                input_videos[inname].append([rmtime, '-', '-', '-', '# remove_part'])
            case ('#OUT', outname, *_):
                pass #outfile = outname
            case ('#DO', cmd, *_): # cmd to run when reading config file
                run(cmd)
            case ('#FINALDO', cmd, *_): # cmd to run after everything is done
                final_cmds.append(cmd)
            case ('#COMMENT', *_):
                pass
            case ('#SET', setting, value, *_):
                args.__setattr__(setting, value)
                print(f'[SET] {setting} = {value}')
            case ('',):
                pass
            case ('#EOF',):
                break
            case _:
                raise ValueError(f'line not understood: \n{line}')

ext = args.tempext
ca, cv = args.ca, args.cv

if cv == 'h264':
    cv += f' -crf {args.crf}'

if outname is None:
    raise ValueError('expected an output name in config. file. example:\n#OUT\toutfile.mp4')

for ini, inname in enumerate(input_videos):
    
    # if len(input_videos) > 9:
    #     raise NotImplementedError('More than 9 clips for one input video: not supported.')

    input_clip = input_videos[inname]
    # dict(zip(['time', 'name', 'cliptime', 'kind', 'comment'], (zip(*input_clip))))
    # time, name, cliptime, kind, *_ = clipinfo
    time, clipnames, cliptimes, kinds, *_ = zip(*input_clip)
    time = [t.split('-') for t in time]
    idx = [[f'{i}-{i+1}', f'{i+1}'] for i in range(len(time))]
    idx = [0] + list(itertools.chain(*idx))
    ts = list(itertools.chain(*time))
    tsec = [parse_time(t) for t in ts]
    assert all([t2 >= t1 for t1, t2 in zip(tsec[:-1], tsec[1:])]), f'time {ts} is not monotonous.'
    
    merge_list = []
    # cut main file
    for idx, tstart, tend in zip(idx, [0] + ts, ts + ['end']):
        if parse_time(tstart) == parse_time(tend) == 0:
            continue
        if parse_time(tstart) == parse_time(tend) == inf:
            continue
        
        start = f'-ss {tstart}' if parse_time(tstart) > 0 else ''
        end = f'-to {tend}' if tend not in ['end', 'END', inf] else ''
        tempname = f'temp/main{ini}.{idx}.{ext}'
        run(f'ffmpeg -accurate_seek -i "{inname}" {start} {end} -c:a {ca} -c:v {cv} -y {tempname}')
        merge_list.append(tempname)
    
    # get clips
    remove_list = []
    for i in range(len(clipnames)):
        clipname = clipnames[i]
        cliptime = cliptimes[i]
        kind = kinds[i]
        clipidx = f'{i}-{i+1}'
        
        if clipname in ['-', '', ' ']:
            remove_list.append(f'temp/clip{ini}.{clipidx}.{ext}')
            continue
                
        if cliptime in ['-', '', ' ']:
            cliptime = ''
        else:
            ts, te = cliptime.split('-')
            cliptime = '%s %s' % (f'-ss {ts}' if ts != '' else '', f'-to {te}' if te != '' else '')
 
        run(f'ffmpeg -accurate_seek -i "{clipname}" {cliptime} -c:a {ca} -c:v {cv} -y temp/clip{ini}.{clipidx}.{ext}')
        if kind == 'a': # replace audio
            raise NotImplementedError
            # print('[WARNING] time of audio must be shoter than the video; not checked')
            copy(f'temp/clip{ini}.{clipidx}.{ext}', f'temp/temp.{ext}')
            run(f'ffmpeg -i temp/temp.{ext} -i temp/main{ini}.{clipidx}.{ext} -c:v h264 -c:a aac ' +
                f'-map 0:a -map 1:v -shortest -y temp/clip{ini}.{clipidx}.{ext}') # 
        elif kind == 'v':
            raise NotImplementedError
            print('[WARNING] time of audio must be shoter than the video; not checked')
            copy(f'temp/clip{ini}.{clipidx}.{ext}', f'temp/temp.{ext}')
            run(f'ffmpeg -i temp/temp.{ext} -i temp/main{ini}.{clipidx}.{ext} -c copy ' +
                f'-map 0:v -map 1:a -y temp/clip{ini}.{clipidx}.{ext}') #  -shortest
        elif kind == 'av':
            pass
        else:
            raise ValueError
        
    # merge videos
    merge_list = [fname.replace('main', 'clip') \
                  if len(fname.split('.')[-2].split('-'))==2 else fname \
                  # if len(fname.split('.')[-2])==2 else fname \
                  for fname in merge_list]
    with open(f'temp/list{ini}.txt', 'w') as f:
        for filename in merge_list:
            if filename in remove_list:
                continue
            filename = filename.replace('temp/', '')
            f.write(f"file '{filename}'\n")
            
    mergename = f'merge{ini}.{ext}'
    run(f'ffmpeg -f concat -i temp/list{ini}.txt -c copy -y temp/{mergename}')
    with open('temp/list.txt', 'w') as f:
        f.write(f"file '{mergename}'\n")

run(f'ffmpeg -f concat -i temp/list.txt -c copy {outname}')
        
for cmd in final_cmds:
    run(cmd)
