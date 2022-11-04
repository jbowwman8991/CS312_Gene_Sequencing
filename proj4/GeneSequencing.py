from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

    def __init__(self):
        self.banded = None
        
    '''
    Align method called by the GUI.
    sequence_1, sequence_2 are the two sequences to be aligned
    banded tells whether or not it should be unrestriced or banded
    align_length is how many base pairs are used in computing the alignment
    '''
    def align(self, sequence_1, sequence_2, banded, align_length):
        self.banded = banded
        self.MaxCharactersToAlign = align_length

        # Adding blank character to the beginning of each sequence
        sequence_1 = ''.join(("-", sequence_1))                                                                             # Time: O(1), Space: O(1)
        sequence_2 = ''.join(("-", sequence_2))                                                                             # Time: O(1), Space: O(1)

        # Creating AlignmentA and AlignmentB
        AlignmentA = ""
        AlignmentB = ""

        # Setting the width and the height
        width = min(align_length + 1, len(sequence_1))                                                                      # Time: O(1), Space: O(1)
        height = min(align_length + 1, len(sequence_2))                                                                     # Time: O(1), Space: O(1)

        # Unrestricted
        if not banded:
            return self.doUnrestricted(sequence_1, sequence_2, align_length, AlignmentA, AlignmentB, width, height)         # Time: O(nm), Space: O(nm)

        # Restricted
        if banded:
            return self.doBanded(sequence_1, sequence_2, align_length, AlignmentA, AlignmentB, width, height)               # Time: O(kn), Space: O(kn)

    def doUnrestricted(self, sequence_1, sequence_2, align_length, AlignmentA, AlignmentB, width, height):                  # Total Time: O(nm), Total Space: O(nm)
        # Initializing main table
        table = [[math.inf for i in range(height)] for j in range(width)]                                                   # Time: O(ij), Space: O(ij)
        # Initializing pointers table that holds pointers to the previous step for backtrace
        pointersTable = [[0 for i in range(height)] for j in range(width)]                                                  # Time: O(ij), Space: O(ij)

        # Setting the first values
        table[0][0] = 0                                                                                                     # Time: O(1), Space: O(1)
        pointersTable[0][0] = 'stop'                                                                                        # Time: O(1), Space: O(1)

        # Setting the base cases
        for i in range(1, width):                                                                                           # Time: O(width), Space: O(1)
            table[i][0] = table[i - 1][0] + INDEL
            pointersTable[i][0] = 'left'                                                                                        
        for j in range(1, height):                                                                                          # Time: O(height), Space: O(1)
            table[0][j] = table[0][j - 1] + INDEL
            pointersTable[0][j] = 'top'

        # Filling the table from the top left
        for i in range(1, width):                                                                                           # Time: O(width), Space: O(1)
            for j in range(1, height):                                                                                      # Time: O(width), Space: O(1)
                left = table[i - 1][j] + INDEL
                top = table[i][j - 1] + INDEL

                # Checking if the letters match
                if sequence_1[i] == sequence_2[j]:                                                                          # Time: O(1), Space: O(1)
                    # Setting value to MATCH
                    diagonal = table[i - 1][j - 1] + MATCH
                else:                                                           
                    # Setting value to sUB
                    diagonal = table[i - 1][j - 1] + SUB
                
                # Getting the minimum
                minimum = min(left, top, diagonal)                                                                          # Time: O(1), Space: O(1)
                table[i][j] = minimum

                # Checking which input was used from left, top, or diagonal
                if minimum == left:                                                                                         # Time: O(1), Space: O(1)
                    pointersTable[i][j] = 'left'
                elif minimum == top:                                                                                        # Time: O(1), Space: O(1)
                    pointersTable[i][j] = 'top'
                elif minimum == diagonal:                                                                                   # Time: O(1), Space: O(1)
                    pointersTable[i][j] = 'diagonal'

        # Creating the alignment cost
        align_cost = table[-1][-1]

        # Getting the alignment
        i = min(len(sequence_1) - 1, align_length)                                                                          # Time: O(1), Space: O(1)
        j = min(len(sequence_2) - 1, align_length)                                                                          # Time: O(1), Space: O(1)

        # While we are not at the beginning 
        while i > 0 or j > 0:                                                                                               # Time: O(ij), Space: O(1)
            # Check what the previous was
            if j > 0 and pointersTable[i][j] == 'top':                                                                      # Time: O(1), Space: O(1)
                # If the previous was 'top' then go up
                AlignmentA = "-" + AlignmentA
                AlignmentB = sequence_2[j] + AlignmentB
                j -= 1
            elif i > 0 and pointersTable[i][j] == 'left':                                                                   # Time: O(1), Space: O(1)
                # If previous was 'left' then go left
                AlignmentA = sequence_1[i] + AlignmentA
                AlignmentB = "-" + AlignmentB
                i -= 1
            elif i > 0 and j > 0 and pointersTable[i][j] == 'diagonal':                                                     # Time: O(1), Space: O(1)
                # If previous was 'diagonal' then go diagonal
                AlignmentA = sequence_1[i] + AlignmentA
                AlignmentB = sequence_2[j] + AlignmentB
                i -= 1
                j -= 1
            else:
                # If none of the above
                return                                                                                                      # Time: O(1), Space: O(1)

        # Returning the cost and the alignments
        return {'align_cost': align_cost, 'seqi_first100': AlignmentA[:100], 'seqj_first100': AlignmentB[:100]}             # Time: O(1), Space: O(1)

    def doBanded(self, sequence_1, sequence_2, align_length, AlignmentA, AlignmentB, width, height):                        # Total Time: O(kn), Total Space: O(kn)
        # Significant length discrepancies can not be aligned
        if abs(width - height) > 100:                                                                                       # Time: O(1), Space: O(1)
            align_cost = math.inf
            alignment_1 = "No Alignment Possible"
            alignment_2 = "No Alignment Possible"
            return {'align_cost': align_cost, 'seqi_first100': alignment_1, 'seqj_first100': alignment_2}                   # Time: O(1), Space: O(1)

        # Resetting the width
        width = 1 + MAXINDELS * 2
        
        # Initializing table
        table = [[float('inf') for i in range(7)] for j in range(height)]                                                   # Time: O(ij), Space: O(ij)
        # Initializing pointers table that holds pointers to the previous step for backtrace
        pointersTable = [[0 for i in range(7)] for j in range(height)]                                                      # Time: O(ij), Space: O(ij)

        # Setting the first values
        table[0][0] = 0                                                                                                     # Time: O(1), Space: O(1)
        pointersTable[0][0] = 'stop'                                                                                        # Time: O(1), Space: O(1)

        # Setting the base cases
        for i in range(1, 4):                                                                                               # Time: O(1), Space: O(1)
            table[0][i] = table[0][i - 1] + INDEL
            pointersTable[0][i] = 'top'

        # Filling the table from the top left
        for i in range(1, height):                                                                                          # Time: O(i), Space: O(1)
            for j in range(width):                                                                                          # Time: O(j), Space: O(1)
                # If values are shifted in the table
                if i > MAXINDELS:                                                                                           # Time: O(1), Space: O(1)
                    # Getting the top value
                    if j + 1 > width - 1:                                                                                   # Time: O(1), Space: O(1)
                        top = math.inf
                    else:                                                               
                        top = table[i - 1][j + 1] + INDEL

                    # Getting the left value
                    if j - 1 < 0:                                                                                           # Time: O(1), Space: O(1)
                        left = math.inf
                    else:
                        left = table[i][j - 1] + INDEL

                    # Adjusting
                    adjust = min(MAXINDELS - i, 0)                                                                          # Time: O(1), Space: O(1)
                    
                    # Checking if the letters match
                    if j - adjust < len(sequence_1) and sequence_1[j - adjust] == sequence_2[i]:                            # Time: O(1), Space: O(1)
                        # Setting value to MATCH
                        diagonal = table[i - 1][j] + MATCH
                    else:
                        # Setting value to SUB
                        diagonal = table[i - 1][j] + SUB

                    # Getting the minimum
                    minimum = min(top, left, diagonal)                                                                      # Time: O(1), Space: O(1)
                    table[i][j] = minimum

                    # Checking which input was used from left, top, or diagonal
                    if minimum == top:                                                                                      # Time: O(1), Space: O(1)
                        pointersTable[i][j] = 'top'
                    elif minimum == left:                                                                                   # Time: O(1), Space: O(1)
                        pointersTable[i][j] = 'left'
                    elif minimum == diagonal:                                                                               # Time: O(1), Space: O(1)
                        pointersTable[i][j] = 'diagonal'

                # If the first three rows don't need an adjustment
                else: 
                    # Getting the top value
                    top = table[i - 1][j] + INDEL
                    # Getting the left value
                    left = table[i][j - 1] + INDEL

                    # Checking if the letters match
                    if sequence_1[i] == sequence_2[j]:                                                                      # Time: O(1), Space: O(1)
                        # Setting value to MATCH
                        diagonal = table[i - 1][j - 1] + MATCH
                    else:
                        # Setting value to SUB
                        diagonal = table[i - 1][j - 1] + SUB

                    # Getting the minimum 
                    minimum = min(top, left, diagonal)                                                                      # Time: O(1), Space: O(1)
                    table[i][j] = minimum

                    # Checking which input was used from left, top, or diagonal
                    if minimum == top:                                                                                      # Time: O(1), Space: O(1)
                        pointersTable[i][j] = 'top'
                    elif minimum == left:                                                                                   # Time: O(1), Space: O(1)
                        pointersTable[i][j] = 'left'
                    elif minimum == diagonal:                                                                               # Time: O(1), Space: O(1)
                        pointersTable[i][j] = 'diagonal'

        # Get the alignment
        i = min(len(sequence_2) - 1, align_length)                                                                          # Time: O(1), Space: O(1)
        j = min(len(sequence_1) - 1, align_length)                                                                          # Time: O(1), Space: O(1)

        # While we are not at the beginning 
        while i > 0 or j > 0:                                                                                               # Time: O(ij), Space: O(1)
            # Checking if the table is shifted 
            if i > MAXINDELS:                                                                                               # Time: O(1), Space: O(1)
                # Adjusting
                adjust = j - i + MAXINDELS
                if i > 0 and j > 0 and pointersTable[i][adjust] == 'left':                                                  # Time: O(1), Space: O(1)
                    # If the previous was 'left' then go left
                    AlignmentA = "-" + AlignmentA
                    AlignmentB = sequence_1[j] + AlignmentB
                    j -= 1
                elif i > 0 and pointersTable[i][adjust] == 'top':                                                           # Time: O(1), Space: O(1)
                    # If the previous was 'top' then go up
                    AlignmentA = sequence_2[i] + AlignmentA
                    AlignmentB = "-" + AlignmentB
                    i -= 1
                elif j > 0 and pointersTable[i][adjust] == 'diagonal':                                                      # Time: O(1), Space: O(1)
                    # If the previous was 'diagonal' then go diagonal
                    AlignmentA = sequence_2[i] + AlignmentA
                    AlignmentB = sequence_1[j] + AlignmentB
                    i -= 1
                    j -= 1
                else:
                    # If none of the above
                    return                                                                                                  # Time: O(1), Space: O(1)
            # If the first three rows don't need an adjustment
            else: 
                if i > 0 and pointersTable[i][j] == 'left':                                                                 # Time: O(1), Space: O(1)
                    # If the previous was 'left' then go left
                    AlignmentA = sequence_2[i] + AlignmentA
                    AlignmentB = "-" + AlignmentB
                    i -= 1
                elif j > 0 and pointersTable[i][j] == 'top':                                                                # Time: O(1), Space: O(1)
                    # If the previous was 'top' then go up
                    AlignmentA = "-" + AlignmentA
                    AlignmentB = sequence_1[j] + AlignmentB
                    j -= 1
                elif i > 0 and j > 0 and pointersTable[i][j] == 'diagonal':                                                 # Time: O(1), Space: O(1)
                    # If the previous was 'diagonal' then go diagonal
                    AlignmentA = sequence_2[i] + AlignmentA
                    AlignmentB = sequence_1[j] + AlignmentB
                    i -= 1
                    j -= 1
                else:
                    # If none of the above
                    return                                                                                                  # Time: O(1), Space: O(1)

        # Swapping AlignmentA with AlignmentB
        AlignmentA, AlignmentB = AlignmentB, AlignmentA                                                                     # Time: O(1), Space: O(1)

        # Getting lenghts
        i = min(len(sequence_2) - 1, align_length)                                                                          # Time: O(1), Space: O(1)
        j = min(len(sequence_1) - 1, align_length)                                                                          # Time: O(1), Space: O(1)

        # Getting the alignment cost
        align_cost = table[i][-abs(j - i) + MAXINDELS]                                                                      # Time: O(1), Space: O(1)

        return {'align_cost': align_cost, 'seqi_first100': AlignmentA[:100], 'seqj_first100': AlignmentB[:100]}             # Time: O(1), Space: O(1)