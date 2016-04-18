# -*- coding: utf-8 -*-
"""
gene_finder parses dna codons.  Currently uses a single string as argument.

@author: Nathan Yee

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq
import time

def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###

complement_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('G')
    'C'
    """
    return complement_dict[nucleotide]

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """

    return ''.join(get_complement(n) for n in dna[::-1])
    new_string = ""

    for i in dna[::-1]: #iterate from the last index to the first index
        new_string = new_string + get_complement(i)

    return new_string

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    #Stop Codons TAA, TAG, TGA

    limit = len(dna) #once we pass limit we stop while loop
    i = 0
    new_string = "" #add dna triples to new_string.  Then return

    while i<limit: #don't loop pase string length

        #first break the loop if we see a stop codon
        if dna[i:i+3] in ['TAA','TAG','TGA']:
            break
        #add dna triple to new_string
        new_string= new_string + dna[i:i+3]

        #increase counter by 3 because we added triple
        i = i + 3

    return new_string


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """

    dna_list = [] #add items to dna_list.  return dna_list at the end
    i = 0 #tells position inside of string dna

    while i < len(dna): #don't go over the length of the list
        temp = ""
        if dna[i:i+3] == 'ATG':
            temp = rest_of_ORF(dna[i:]) #temporary variable append and add length to i
            dna_list.append(temp)
        i = i + len(temp) + 3

    return dna_list


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """

    dna_list = []
    for i in range(3): #go through each possible frames
        dna_list = dna_list + find_all_ORFs_oneframe(dna[i:])
    return dna_list


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """

    string_list = [dna, get_reverse_complement(dna)] #string_list contains both dna strands
    dna_list = [] #add orf's to list then return


    for list in string_list:
        temp_list = find_all_ORFs(list)
        for item in temp_list:
            dna_list.append(item)

    return dna_list

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    strands = find_all_ORFs_both_strands(dna) #get the strands for a string of dna
    longest_strand = "" #declare longest_strand variable

    for i in strands: #go through all the strands
        if len(longest_strand) < len(i): #compares length of current longest strand do index
            longest_strand = i #assign new longest strand if index is greater than previous
    return longest_strand #return the longest strand as a string

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    maximum_length = 0 #maximum length will store the strand with the maximum length
    string = "" #this is the string that we will return
    for i in range(num_trials): #run program num_trials times
        shuffled = longest_ORF(shuffle_string(dna)) #get the longest string frum shuffled input
        if len(shuffled) > maximum_length: #compare length of string to previous maximum
            maximum_length = len(shuffled) #assign new maximum if true
            string = shuffled
    return string

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    protein = "" #protein stores each amino value for triples
    for i in range(0, len(dna)-2, 3):
        protein.join(aa_table[dna[i:i+3]])
    return protein

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    AA_list = []
    threshold = longest_ORF_noncoding(dna,1500)
    all_ORFs = find_all_ORFs_both_strands(dna)
    for item in all_ORFs:
        if len(item) > len(threshold):
            AA_list.append(coding_strand_to_AA(item))
    return AA_list

start_time = time.time()
from load import load_contigs 
contigs = load_contigs()
name,dna = contigs[5]
dna = dna[1:10000]

gene_finder(dna)
print("---%s seconds ---" % (time.time()-start_time))

if __name__ == "__main__":
    import doctest
    doctest.testmod()
