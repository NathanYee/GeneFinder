# -*- coding: utf-8 -*-
"""
gene_finder parses dna codons.  Currently uses a single string as argument.

@author: Nathan Yee

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


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

    #series of if elif statements that return complementary nucleotides
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    else:
        print "Invalid Arguemnt not a nucleotide"
        return None

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
    for i in range(num_trials): #run program num_trials times
        shuffled = longest_ORF(shuffle_string(dna)) #get the longest string frum shuffled input
        if len(shuffled) > maximum_length: #compare length of string to previous maximum
            maximum_length = len(shuffled) #assign new maximum if true
    return maximum_length #return the length as an int

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
    i = 0
    protein = "" #protein stores each amino value for triples
    amino = "" #a value which we append to protein, will change depending on triple catigorization
    for i in range(0,len(dna),3): #we step 3 because we want triples
        triple = dna[i:i+3] #get the next three codons from dna

        if len (triple) == 3: #make sure triple has length 3, appending occures at the end of this statement
            #lots of if elif statements to categorize dna to amino acid
            #each will assign a new value to amino
            if triple in ["TTT","TTC"]:
                amino = "F"
            elif triple in ['TTA','TTG','CTT','CTC','CTA','CTG']:
                amino = "L"
            elif triple in ['ATT','ATC','ATA']:
                amino = "I"
            elif triple in ["ATG"]:
                amino = "M"
            elif triple in ['GTT','GTC','GTA','GTG']:
                amino = "V"
            elif triple in ['TCT','TCC','TCA','TCG']:
                amino = "S"
            elif triple in ['CCT','CCC','CCA','CCG']:
                amino = "P"
            elif triple in ['ACT','ACC','ACA','ACG']:
                amino = "T"
            elif triple in ['GCT','GCC','GCA','GCG']:
                amino = "A"
            elif triple in ['TAT','TAC']:
                amino = "Y"
            elif triple in ['CAT','CAC']:
                amino = "H"
            elif triple in ['CAA','CAG']:
                amino = "Q"
            elif triple in ['AAT','AAC']:
                amino = "N"
            elif triple in ['AAA','AAG']:
                amino = "K"
            elif triple in ['GAT','GAC']:
                amino = "D"
            elif triple in ['GAA','GAG']:
                amino = "E"
            elif triple in ['TGT','TGC']:
                amino = "C"
            elif triple in ['TGG']:
                amino = "W"
            elif triple in ['CGT','CGC','CGA','CGG']:
                amino = "R"
            elif triple in ['AGT','AGC']:
                amino = "S"
            elif triple in ['AGT','AGC']:
                amino = "R"
            elif triple in ['GGT','GGC','GGA','GGG']:
                amino = "G"
            # end categorization
            protein = protein + amino #append amino value to protein
    return protein #return string

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """

if __name__ == "__main__":
    import doctest
    doctest.testmod()
