# Hangman game
#

# -----------------------------------
# Helper code
# You don't need to understand this helper code,
# but you will have to know how to use the functions
# (so be sure to read the docstrings!)

import random

WORDLIST_FILENAME = "words.txt"

def loadWords():
    """
    Returns a list of valid words. Words are strings of lowercase letters.
    
    Depending on the size of the word list, this function may
    take a while to finish.
    """
    # print("Loading word list from file...")
    # inFile: file
    inFile = open(WORDLIST_FILENAME, 'r')
    # line: string
    line = inFile.readline()
    # wordlist: list of strings
    wordlist = line.split()
    # print("  ", len(wordlist), "words loaded.")
    return wordlist

def chooseWord(wordlist):
    """
    wordlist (list): list of words (strings)

    Returns a word from wordlist at random
    """
    return random.choice(wordlist)

# end of helper code
# -----------------------------------


def isWordGuessed(secretWord, lettersGuessed):
    '''
    secretWord: string, the word the user is guessing
    lettersGuessed: list, what letters have been guessed so far
    returns: boolean, True if all the letters of secretWord are in lettersGuessed;
      False otherwise
    '''
    # FILL IN YOUR CODE HERE...
    guessedFlag = True
    # Convert string into list
    secretWordChars = list(secretWord)
    for char in secretWordChars    :
        if char not in lettersGuessed:
            guessedFlag = False
    return guessedFlag 


def getGuessedWord(secretWord, lettersGuessed):
    '''
    secretWord: string, the word the user is guessing
    lettersGuessed: list, what letters have been guessed so far
    returns: string, comprised of letters and underscores that represents
      what letters in secretWord have been guessed so far.
    '''
    # FILL IN YOUR CODE HERE...
    secretWordChars = list(secretWord)
    outputChars = ''
    for char in secretWordChars    :
        if char not in lettersGuessed:
             outputChars = outputChars + '_ '
        else: 
            outputChars = outputChars + char
    return outputChars   


def getAvailableLetters(lettersGuessed):
    '''
    lettersGuessed: list, what letters have been guessed so far
    returns: string, comprised of letters that represents what letters have not
      yet been guessed.
    '''
    # FILL IN YOUR CODE HERE...
    alphabeth = list('abcdefghijklmnopqrstuvwxyz')
    available = ''
    for char in alphabeth:
        if char not in lettersGuessed:
            available = available + char
    return available      

def hangman(secretWord):
    '''
    secretWord: string, the secret word to guess.

    Starts up an interactive game of Hangman.

    * At the start of the game, let the user know how many 
      letters the secretWord contains.

    * Ask the user to supply one guess (i.e. letter) per round.

    * The user should receive feedback immediately after each guess 
      about whether their guess appears in the computers word.

    * After each round, you should also display to the user the 
      partially guessed word so far, as well as letters that the 
      user has not yet guessed.

    Follows the other limitations detailed in the problem write-up.
    '''
    print("")
    print("------------")
    print("Welcome to the game Hangman!")
    print("------------")
    print("")
    print("I am thinking of a word that is ",  len(secretWord)," letters long.")   
    maxNGuesses = 8
    guessedFlag = False
    nGuesses = 0
    lettersGuessed = ''
    while not(guessedFlag) and nGuesses < maxNGuesses:
        print("")
        print("You have ", maxNGuesses - nGuesses, " guesses left.")
        availableLetters = getAvailableLetters(lettersGuessed)
        print("Available letters:  "+availableLetters)
        guessCheck = False
        # Ask for a guess, check if letter guessed is OK (not in the ones already guessed)
        while not(guessCheck):            
            letterGuess = input("Please guess a letter: ")
            guessInLowerCase = letterGuess.lower()
            if guessInLowerCase in lettersGuessed:
                guessCheck = False
                print("Oops! You've already guessed that letter: ", getGuessedWord(secretWord, lettersGuessed))
                print("------------")
            else: 
                lettersGuessed = lettersGuessed+letterGuess 
                guessCheck = True
        # Check if the provided letter is part of the word to guess.
        # Print the word as guessed so far
        # If the guess is wrong, increment the number of guesses       
        secretWordChars = list(secretWord)
        if letterGuess in secretWordChars:
            print("Good guess: ", getGuessedWord(secretWord, lettersGuessed))
        else:
            print("Oops! That letter is not in my word: ", getGuessedWord(secretWord, lettersGuessed))
            nGuesses += 1
        print("------------")
        # Check if the word has been guessed
        # If it is, print message of congrats
        guessedFlag = isWordGuessed(secretWord, lettersGuessed)
        if guessedFlag:
            print("Congratulations, you won!")
        if nGuesses == maxNGuesses:
            print("Sorry, you ran out of guesses. The word was "+secretWord+"." )
                  
# Play hangman!!
# Load the list of words into the variable wordlist
# so that it can be accessed from anywhere in the program
wordlist = loadWords()
secretWord = chooseWord(wordlist)
#print(secretWord)
hangman(secretWord)
