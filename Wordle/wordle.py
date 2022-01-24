# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 15:37:15 2022

@author: Yuchen Wang
"""

from string import ascii_lowercase
import numpy as np

doc = '''\
This is a simple script that can help you play the Wordle game.

There are 3 modes:
Help mode (h): Got stuck? Use this mode to input all your tries (and the colors) so far, and I will suggest a word for you.
Computer mode (c): The computer will try to guess the word. All you need to input is the color of the computer\'s try.
Play (p): You can play a Wordle game here.

Notes on colors:
Use g, y, b for green, yellow, black (grey). Example:gbybb.
Green (g): This letter is in the word and in the correct spot.
Yellow (y): This letter is in this word but in the wrong spot.
Black (b): This letter is not in this word in any spot.
'''

ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])

with open('engword.txt') as f:
    words = eval(f.read())
wordlist = np.array(words)
words = [[i for i in word] for word in words]
words = np.array(words)
words_idx = words.copy().astype('<U2')

for i, c in enumerate(ascii_lowercase):
    words_idx[np.where(words_idx==c)] = i
words_idx = words_idx.astype(int)

cs = np.array([i for i in ascii_lowercase])
letter_occur = [np.isin(cs, word) for word in words]
letter_occur = np.stack(letter_occur)
let_num = np.sum(letter_occur, axis=0)

let_unique = [np.unique(idx) for idx in words_idx]
word_score = np.array([np.sum(let_num[word]) for word in let_unique])

print(doc)
# print('Note: use g, y, b for green, yellow, black (grey).\nexample:gbybb')

while True:
    mode = input('"h" for help mode, "c" for computer mode, "p" for play, "e" to exit. >>> ')
    if mode in ['e', 'exit', 'q', 'quit']:
        break
    guess_n = 6
    # guesses = []
    if mode == 'p':
        word = np.random.choice(wordlist)
        for i in range(guess_n):
            while True:
                guess = input(f'Your {ordinal(i+1)} try: >>> ')
                if guess in wordlist:
                    break
                else:
                    print(f'Word "{guess}" does not exit or is not a five letter word.')
            color = ['-', '-', '-', '-', '-']
            lets = []
            if guess == word:
                print(f'Congratulations! The answer is "{word}".')
                break
            elif guess != word and i+1 == guess_n:
                print(f'You lose! The word is {word}.')
                break
            for j, c in enumerate(guess):
                if c in word and c not in lets:
                    if word[j] == c:
                        color[j] = 'g'
                    else:
                        color[j] = 'y'
                    lets.append(c)
            for j, c in enumerate(guess):
                if c not in word and c not in lets:
                    color[j] = 'b'
            print(f'Your color is {"".join(color)}.')
    elif mode in ['c', 'h']:
        word_left = np.full(len(words), True)
        for i in range(guess_n):
            onlyone = 'This is the only word possible!!' if np.sum(word_left) == 1 else ''
            if mode == 'c':
                guess = wordlist[word_left][np.argmax(word_score[word_left])]
                print(f'The {ordinal(i+1)} guess is "{guess}". {onlyone}')
            elif mode == 'h':
                while True:
                    guess = input(f'input your {ordinal(i+1)} guess, or "s" to suggest one for you: >>> ')
                    if guess == 's':
                        print('My suggestion is:', wordlist[word_left][np.argmax(word_score[word_left])]+'.', onlyone)
                    elif guess in wordlist:
                        break
                    else:
                        print(f'Word "{guess}" does not exit or is not a five letter word.')
            else:
                raise ValueError(f'Unknown mode "{mode}".')
            # guesses.append(guess)
            while True:
                color = input('Input color (e.g. gbybb) for this guess: >>> ')
                if len(color) == 5 and all([co in ['b', 'y', 'g'] for co in color]):
                    break
                else:
                    print('Color should be a string with 5 letters in (\'b\', \'y\', \'g\').')
            if color == 'ggggg':
                print('Hooray!')
                break
            else:
                lets = []
                for j, col in enumerate(color):
                    if col == 'g':
                        lets.append(guess[j])
                        word_left &= words[:, j] == guess[j]
                    elif col == 'y':
                        lets.append(guess[j])
                        word_left &= words[:, j] != guess[j]
                        word_left &= np.array([guess[j] in word for word in wordlist])
                for j, col in enumerate(color):
                    if col == 'b' and guess[j] not in lets:
                        word_left &= np.array([guess[j] not in word for word in wordlist])
                if np.sum(word_left) == 0:
                    print('No possible found! Check if your input is correct.')
                    break        
    else:
        print(f'Unknown mode "{mode}".')
