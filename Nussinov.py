import numpy as np
import RNA

complementaire={"C":"G","G":"C","A":"U","U":"A" }

def Appariement(structure: str, L=None)->list:
    """genère l'appariement d'une structure à partir de sa forme '.((...-()..))

    Args:
        structure (str): la structure
        L ([type], optional): La liste des appariements. Defaults to None.

    Returns:
        list[int]: liste telle que L[i]=j si i est apparayé à j ou -1 sinon
    """
    L = [-1]*len(structure)
    stack = []
    for i in range(len(structure)):
        if structure[i]=='(':
            stack.append(i)
        elif structure[i]==')':
            if len(stack)==0:
                raise Exception(") en trop")
            j = stack.pop()
            L[i] = j
            L[j] = i
    if len(stack)>0:
        raise Exception("( en trop")
    return L

def Structure(base_pair: list, l: int):
    """Generate the tructure ((..().)) from a base pair list 

    Args:
        base_pair (list[(int, int)]): base pair list
        l (int): structure lenght

    Returns:
        str: structure
    """

    stru=["."]*l
    for coord in base_pair:
        stru[coord[0]-1]="("
        stru[coord[1]-1]=")"
    return ("".join(stru))

class Simple_Nussinov:
    """Class embeding all the function to compute Nussinov algorithm on one sequence
    """

    def __init__(self, sequence: str, m=0):
        """Constructor of the class

        Args:
            sequence (str): sequence string
            m (int, optional): the minimum lenght of folds. Defaults to 0.
        """
        self.sequence = sequence.replace('-', '')
        self.m = m
        n = len(self.sequence)
        self.M = np.zeros([n+1, n+1])
        self.base_pair = []
        self.structure = ''

    def Compute(self, print_loading=False):
        """Compute the Nussinov algorithm to fill the matrix M

        Returns:
            np.array: the matrix
        """

        i,j=1,2
        c = 2
        n = len(self.sequence)
        count = 0
        count_max = n*(n-1)//2
        while (c<n+1):
            count+=1
            if print_loading:
                print("="*int(count/count_max*10) + "-"*int((1-count/count_max))*10, end='\r')
            pair=[]
            for k in range(i, j-self.m):
                if(k<0):
                    break
                else :
                    if(complementaire[self.sequence[j-1]] == self.sequence[k-1]):
                        pair.append(self.M[i][k-1]+1+self.M[k+1][j-1])
                if(pair==[]):
                    p=0
                else : 
                    p=max(pair)
                self.M[i][j]=max(self.M[i][j-1], p)

            if (j==n):
                i=1
                c+=1
                j=c
            else :
                i,j = i+1, j+1
        if print_loading:
            print("Done                ")
        return self.M

    def __traceback(self, i=None, j=None):
        """Reccursive traceback method

        Args:
            i (int, optional): i index. Defaults to None.
            j (int, optional): j index. Defaults to None.

        Returns:
            list[(int, int)]: base pair
        """

        if i==None:
            i=1
            j=len(self.M[0])-1
            self.base_pair = []

        if(j<=i):
            return self.base_pair
        elif(self.M[i][j]==self.M[i][j-1]):
            self.__traceback(i, j-1)
            return self.base_pair
        else:
            for k in range(i,j):
                if(complementaire[self.sequence[j-1]] == self.sequence[k-1]):
                    if(self.M[i][j] == self.M[i][k-1] + self.M[k+1][j-1] +1):
                        self.base_pair.append((k,j))
                        self.__traceback(i, k-1)
                        self.__traceback(k+1, j-1)
                        return self.base_pair

    def Compute_Traceback(self):
        """Compute the traceback method
        """
        self.__traceback()
        self.structure = Structure(self.base_pair, len(self.sequence))
        
class ACC_Computation:
    """Class embeding the methods to compute the ACC vectors
    """

    def __init__(self, structures, appariements, method='MTA', alignments=None, lbda=1):
        """Constructor of the class

        Args:
            structures (list[str]): structures of the sequences ((..()))
            appariements (list[(int, int)]): appariements of the sequences
            method (str, optional): MTA or MEA. Defaults to 'MTA'.
            alignments (list[str], optional): structures with gaps. Defaults to None.
            lbda (int, optional): lambda to compute PAIRED_ACC. Defaults to 1.
        """
        
        self.structures = structures
        self.appariements = appariements
        self.alignments = alignments
        self.method = method
        self.lbda = lbda
        self.UNPAIRED_ACC = np.zeros(len(structures[0]))
        self.PAIRED_ACC = np.zeros([len(structures[0]), len(structures[0])])
        self.bpp_list = []
        self.unp_list = []

        if self.method == 'MEA':
            md = RNA.md()
            # activate unique multibranch loop decomposition
            md.uniq_ML = 1
            # create fold compound object
            fc_list = [RNA.fold_compound(sequence, md) for sequence in self.alignments]
            # compute MFE
            
            for i, fc in enumerate(fc_list):
                (ss, mfe) = fc.mfe()
                # rescale Boltzmann factors according to MFE; rescaling avoids numerical problems for long sequences
                fc.exp_params_rescale(mfe)
                # compute partition function to fill DP matrices
                fc.pf()

                # get a matrix of base pair probabilities (attention: this matrix is 1-based), only entries i<j are meaningful
                bpp = fc.bpp()
                self.bpp_list.append(bpp)
                unp = [1-sum(bpp[i+1]) for i in range(0,len(self.alignments[i]))]
                self.unp_list.append(unp)

    def Compute_UNPAIRED_ACC(self):
        """Compute the UNPAIRED_ACC vector

        Returns:
            list[int]: UNPAIRED_ACC vector
        """

        if self.method == 'MTA':
            K, N = len(self.structures), len(self.structures[0])
            for j in range(N):
                s = 0
                for k in range(K):
                    if self.method=='MTA' and self.structures[k][j]=='.': s+=1
                    elif self.method=='MEA': s+= self.unp_list[k][j+1]
                self.UNPAIRED_ACC[j]=s
            
        return self.UNPAIRED_ACC

    def Compute_PAIRED_ACC(self):
        """Compute PAIRED_ACC matrix

        Returns:
            np.array: PAIRED_ACC matrix
        """
        K, N = len(self.structures), len(self.structures[0])
        for i in range(N):
            for j in range(N):
                s=0
                for k in range(K):
                    if self.method=='MTA' and self.appariements[k][i]==j:s+=1
                    elif self.method=='MEA': s+= self.bpp_list[k][i+1][j+1]
                self.PAIRED_ACC[i][j] = 2*s*self.lbda
        return self.PAIRED_ACC

    def Add_Zeros(self):
        """Add a zero line and column to PAIRED_ACC matrix and a 0 to UNPAIRED_ACC
        """

        self.PAIRED_ACC = np.insert(self.PAIRED_ACC,0, [0]*self.PAIRED_ACC.shape[0], axis = 0)
        self.PAIRED_ACC = np.insert(self.PAIRED_ACC,0, [0]*self.PAIRED_ACC.shape[0], axis = 1)
        self.UNPAIRED_ACC = np.insert(self.UNPAIRED_ACC, 0, [0], axis=0)

    def Remove_Zeros(self):
        """Opposite as Add_Zeros
        """

        self.PAIRED_ACC = np.delete(self.PAIRED_ACC, 0, 0)
        self.PAIRED_ACC = np.delete(self.PAIRED_ACC, 0, 1)
        self.UNPAIRED_ACC = np.delete(self.UNPAIRED_ACC, 0, 0)

class Multiple_Nussinov:
    """Class embeding all the method to compute the Nussinov variant to find a consensus
    """
    def __init__(self, alignments: list, structures: list, m: int, method = 'MTA') -> None:
        """Constructor of the class

        Args:
            alignments (list[str]): structures with gaps
            structures (list[str]): structures ((..()))
            m (int): minimum fold lenght
            method (str, optional): MEA or MTA. Defaults to 'MTA'.
        """
        
        self.alignments = alignments
        self.structures = structures
        self.m = m
        self.method = method
        K, N = len(self.alignments), len(self.alignments[0])
        self.M = np.zeros((N,N))
        for i in range(K):self.structures[i] = self.structures[i].replace('-','')
        for n in range(K):
            for i in range(N):
                if self.alignments[n][i]=='-':
                    self.structures[n] = self.structures[n][0:i] + '-' + self.structures[n][i:]
        self.appariements = [Appariement(s) for s in self.structures]
        self.ACC = ACC_Computation(self.structures, self.appariements, method=self.method, alignments=self.alignments)
        self.base_pair = []
        self.structure = ''

    def ACC_Compute(self):
        """Compute ACC
        """

        self.ACC.Compute_PAIRED_ACC()
        self.ACC.Compute_UNPAIRED_ACC()
    
    def Compute_Nussinov_ACC(self, print_loading=False)->np.array:
        """Compute the Nussinov variant to calculate the MTA

        Returns:
            np.array: The nussinov matrix
        """
        K, N = len(self.alignments), len(self.alignments[0])


        UNPAIRED_ACC_1 = np.zeros(N+1)
        for i in range(N):
            UNPAIRED_ACC_1[i+1] = self.ACC.UNPAIRED_ACC[i]
        PAIRED_ACC_1 = np.zeros([N+1, N+1])
        for i in range(N):
            for j in range(N):
                PAIRED_ACC_1[i+1, j+1] = self.ACC.PAIRED_ACC[i,j]
        i,j=0,0
        c = 0
        count = 0
        count_max = (N-1)*N//2
        while (c<N):
            if print_loading:
                print("="*int(count/count_max*10) + "-"*int((1-count/count_max)*10), end='\r')
            count+=1
            pair=[self.M[i][j-1] + self.ACC.UNPAIRED_ACC[j]]
            for k in range(i, j-self.m):
                if(k<0):
                    break
                else :
                    p = self.ACC.PAIRED_ACC[k][j]
                    if p > 0:
                        #print(str(k) + " lié à " + str(j))
                        pair.append(self.M[i][k-1]+p+self.M[k+1][j-1])
            p=max(pair)
            self.M[i][j]=p

            if (j==N-1):
                i=0
                c+=1
                j=c
            else :
                i,j = i+1, j+1
        if print_loading:
            print("Done          ")  
        
        return(self.M)

    def __traceback_sequence(self, i=None, j=None):
        """Reccursive traceback method

        Args:
            i (int, optional): i index. Defaults to None.
            j (int, optional): j index. Defaults to None.

        Returns:
            list[(int, int)]: base pair
        """
        if i==None:
            i=1
            j=len(self.M[0])-1
            self.base_pair = []

        if(j<=i):
            return self.base_pair
        elif(self.M[i][j]==self.M[i][j-1] + self.ACC.UNPAIRED_ACC[j]):
            self.__traceback_sequence(i, j-1)
            return self.base_pair
        else:
            for k in range(i,j):
                p = self.ACC.PAIRED_ACC[k][j]
                if p>0:
                    if(self.M[i][j] == self.M[i][k-1] + self.M[k+1][j-1] +p):
                        self.base_pair.append((k,j))
                        self.__traceback_sequence(i, k-1)
                        self.__traceback_sequence(k+1, j-1)
                        return self.base_pair
    
    def Compute_Traceback(self):
        """Compute the traceback algorithm
        """

        self.M = np.insert(self.M, 0, [0]*self.M.shape[0], axis=0)
        self.M = np.insert(self.M, 0, [0]*self.M.shape[0], axis=1)
        self.ACC.Add_Zeros()
        self.__traceback_sequence()
        self.M = np.delete(self.M, 0, 1)
        self.M = np.delete(self.M, 0, 0)
        self.ACC.Remove_Zeros()
        self.structure = Structure(self.base_pair, len(self.alignments[0]))
