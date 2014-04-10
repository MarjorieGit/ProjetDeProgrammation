# Code recupere : https://github.com/alevchuk/pairwise-alignment-in-python/blob/master/alignment.py


def zeros(form):
    matrix = []
    for x in range(form[0]):
        matrix.append([])
        for y in range(form[1]):
            matrix[-1].append(0)
    return matrix

match_award = 10 # quand les paires sont complementaires AT CG
mismatch_penalty = -5 # quand les paires de nt sont identiques
mismatch_adenin = -4 # quand les paires sont differentes et associees A-C-G
mismatch_uridin = -3 # quand les paires sont differentes et associees T-C-G
gap_penalty = -7 # both for opening and extanding


def match_score(alpha, beta):
    if alpha == beta:
        return mismatch_penalty
    elif alpha == '-' or beta == '-':
        return gap_penalty
    elif (alpha == 'A' and beta =='C') or (alpha == 'A' and beta =='G') or \
    (alpha == 'C' and beta =='A') or (alpha == 'G' and beta =='A'):
        return mismatch_adenin
    elif (alpha == 'T' and beta =='C') or (alpha == 'T' and beta =='G') or \
    (alpha == 'C' and beta =='T') or (alpha == 'G' and beta =='T'):
        return mismatch_uridin
    else:
        return match_award

def finalize(align1, align2):
    align1 = align1[::-1] #reverse sequence 1
    align2 = align2[::-1] #reverse sequence 2
    
    i,j = 0,0
    
    #calcuate identity, score and aligned sequeces
    symbol = ''
    found = 0
    score = 0
    #identity = 0
    for i in range(0,len(align1)):
        # if two AAs are the same, then output the letter
        if align1[i] == align2[i]:
            symbol += ' '
            #identity = identity + 1
            score += match_score(align1[i], align2[i])
    
        elif (align1[i] == 'A' and align2[i] == 'T') or (align1[i] == 'G' and align2[i] == 'C') or (align1[i] == 'T' and align2[i] == 'A') or (align1[i] == 'C' and align2[i] == 'G'):
            score += match_score(align1[i], align2[i])
            symbol += '|'
            # if they are not identical and none of them is gap
    
        elif align1[i] != align2[i] and align1[i] != '-' and align2[i] != '-':
                score += match_score(align1[i], align2[i])
                symbol += ' '
                found = 0
            
        
    
        #if one of them is a gap, output a space
        elif align1[i] == '-' or align2[i] == '-':
            symbol += ' '
            score += gap_penalty
    
    #identity = float(identity) / len(align1) * 100
    
    #print 'Identity =', '%3.3f' % identity, 'percent'
    print 'Score =', score
    print align1
    print symbol
    print align2
    dict = {'align1' : align1, 'symbol': symbol, 'align2' : align2}
    return dict

def Transformation_TU(seq) :
    seq = seq.replace('U','T')
    return seq


def Smith_Waterman(seq1, seq2):
    m, n = len(seq1), len(seq2) # length of two sequences
    
    # Generate DP table and traceback path pointer matrix
    score = zeros((m+1, n+1)) # the DP table    DP = programmation dynamique, transformation seq1 en seq2
    pointer = zeros((m+1, n+1)) # to store the traceback path
    
    max_score = 0 # initial maximum score in DP table
    # Calculate DP table and mark pointers
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            score_diagonal = score[i-1][j-1] + match_score(seq1[i-1], seq2[j-1])
            score_up = score[i][j-1] + gap_penalty
            score_left = score[i-1][j] + gap_penalty
            score[i][j] = max(0,score_left, score_up, score_diagonal)
            if score[i][j] == 0:
                pointer[i][j] = 0 # 0 means end of the path
            if score[i][j] == score_left:
                pointer[i][j] = 1 # 1 means trace up
            if score[i][j] == score_up:
                pointer[i][j] = 2 # 2 means trace left
            if score[i][j] == score_diagonal:
                pointer[i][j] = 3 # 3 means trace diagonal
            if score[i][j] >= max_score:
                max_i = i
                max_j = j
                max_score = score[i][j];
    
    align1, align2 = '', '' # initial sequences
    
    i,j = max_i,max_j # indices of path starting point
    
    #traceback, follow pointers
    while pointer[i][j] != 0:
        if pointer[i][j] == 3:
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1
            j -= 1
        elif pointer[i][j] == 2:
            align1 += '-'
            align2 += seq2[j-1]
            j -= 1
        elif pointer[i][j] == 1:
            align1 += seq1[i-1]
            align2 += '-'
            i -= 1

    dict = finalize(align1, align2)
    return dict

mRNA = 'CTTTTCCACATGCGTGCTGCGTGGATATGTCGCAATCCCGGTGCTATGATCCCTGACCTACGGCGCCTCGCGCCCCGCTCCAAGCAGCACCCGAACCGGCCGCGCCCATGTCGCGTGTCCTACGACCTGATGAAAATGGAGAGTTGAATCATGGCTTTTGAACTGCCGCCGCTGCCGTACGAGAAAAACGCCCTCGAGCCG'
# reverse du sRNA
sRNA = 'ttgaccaacaaatgagaagcgttattatttgaaaactgatcgcgagatcagcggtaaactgacagagcggctactcttcaggttatctcctcatcaggctaattacggttttcgacccggctccgccgggtctttttttg'
sRNA = sRNA[::-1] # reverse
mRNA = mRNA.upper()
sRNA = sRNA.upper()

# traitement des sequences pour transormer les U en T si besoin
mRNA = Transformation_TU(mRNA)
sRNA = Transformation_TU(sRNA)
# d'abord mRNA puis sRNA    
Smith_Waterman(mRNA,sRNA)
