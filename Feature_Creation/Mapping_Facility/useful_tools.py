#Author: Kerstin Lenhof
#Date: 6.11.2018
#Useful tools for my project ...


#Case insensitive comparison normalizer 
import unicodedata

def normalize_caseless(word):
	#return word.casefold()
	return unicodedata.normalize("NFKD", word.casefold())

def equal_caseless(word1, word2):
	return normalize_caseless(word1) == normalize_caseless(word2)