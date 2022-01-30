# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2


    What this function is doing:
    1. load each species BRD2 fasta file
    2. Aligning each species BRD2 to human BRD2. Saving the results of the alignment in a variable. The results of the alignment are in the form of a tuple
    of the form (alignment score, human sequence alignment, other species sequence alignment)
    3. Appending the name of the species to the end of the tuple
    4. Adding each alignment result to a list, `alignment_list`
    5. Sort this list of tuples -- where each tuple contains the alignment score, human seq alignment, other species seq alignment, and species name (in that order) -- by
    sorting based on the 0th element in each tuple, the alignment score, in descending order. Now the list contains the tuples in order of most to least similar to humans because the alignment
    scores are in descending order
    6. create a new list of tuples, each with only the species name and alignment score for each species, ordered from most to least similar to humans
    7. Print this list of tuples
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")


    alignment = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat",-10,-1)

    #NeedlemanWunsch.align returns a tuple of form (alignment score, seqA alignment, seqB alignment)
    gg_alignment = alignment.align(hs_seq,gg_seq) + ('Gallus_gallus',) #Add species name to end of alignment tuple to later print species names in order of similarity
    mm_alignment = alignment.align(hs_seq,mm_seq) + ('Mus_musculus',)
    br_alignment = alignment.align(hs_seq,br_seq) + ('Balaeniceps_rex',)
    tt_alignment = alignment.align(hs_seq,tt_seq) + ('tursiops_truncatus',)

    def sort_key(list_alignments):
        return list_alignments[0] #create a key to sort each alignment tuple by the alignment score, which is the 0th index

    alignment_list = [gg_alignment,mm_alignment,br_alignment,tt_alignment]
    alignment_list.sort(key = sort_key,reverse=True) #sort the list  of alignment tuples in descending order by using the alignment score (0th index) in each tuple. This sorts in place
    species_similarity = [(tup[-1],tup[0]) for tup in alignment_list] #take alignment list that has been sorted by alignment score and put the species names and similarity score for each species in a new list
    


    #print each species, plus their alignment score, in order of most to least similar to human
    print(f"From most similar to least similar, the species most similar to human, plus the alignment score, is:\n{species_similarity}")
   

if __name__ == "__main__":
    main()
