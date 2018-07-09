import numpy as np
import pandas as pd
from scipy.spatial import distance
from typing import Any, Dict, List, Tuple
from ete3 import Tree

class UPGMA:
    def newickFromListOfLists(lists_allthewaydown: List[Any]) -> str:
        """ 
            lists_allthewaydown: a list of lists (possibly of more lists)
            returns: the representation of the list of lists in newick form
        """
        return repr(lists_allthewaydown).replace(' ','').replace("'",'').replace('[','(').replace(']',')')

    
    def clusterFromListOfLists(lists_allthewaydown: List[Any]) -> List[str]:
        """ 
            lists_allthewaydown: a list of lists (possibly of more lists)
            returns: a list of clusters
        """
 
        return [''.join(c for c in repr(s) if c.isalpha()) for s in lists_allthewaydown]

    def computeDifference(traits: Any, i: int, j: int) -> int:
        """
            traits: a trait matrix
            returns: the Hamming distance (number of differences) between the ith and jth columns
        """
        return traits.shape[0]*distance.hamming(traits[:,i], traits[:,j])
    
    def diffMatrixBuilder(traits: Any) -> Any:
        """
            traits: a trait matrix
            returns: a difference matrix based on the trait matrix
        """
        diffs = np.zeros((traits.shape[0], traits.shape[0]), dtype=float)
        for i in range(diffs.shape[0]):
            for j in range(diffs.shape[0]):
                diffs[i,j] = UPGMA.computeDifference(traits, i, j)
            diffs[i,i] = np.inf
        return diffs
    
    def outputMatrix(output_file: str, diffs: Any, clusters: List[Any]) -> None:
        """
            output_file: A file to which the new matrix will be appended
            diffs: A differences np matrix that needs to be displayed
            clusters: A list of clusters thusfar. Length should equal the
                number of columns/rows in the differences matrix. Each cluster
                is a list of labels
        """

        out = open(output_file, 'a')
        names = UPGMA.clusterFromListOfLists(clusters)
        df = pd.DataFrame(diffs, index=names, columns=names)
        pd.options.display.float_format = '{:,.2f}'.format
        out.write(df.to_string()+'\n'*2)
        print("Intermediate Matrix: \n", df, '\n')
        
    def upgmaStep(diffs: Any, clusters: List[Any]) -> Any:
        """
            diffs: An mxm numpy array of distances between the m clusters at this step
            clusters: A length-m list consisting of the list of elements in each cluster at this step
            returns: An (m-1)x(m-1) array of distances between the (m-1) clusters after
                collapsing another cluster
        """
        """
            Given a difference matrix, takes the minimum value of this matrix
            and uses it to compute and output the next difference matrix
        """

        i,j = np.unravel_index(np.argmin(diffs, axis=None), diffs.shape)
        s_i = len(UPGMA.clusterFromListOfLists(clusters[i]))
        s_j = len(UPGMA.clusterFromListOfLists(clusters[j]))

        new_column = (s_i*diffs[:,i] + s_j*diffs[:,j])/(s_i+s_j)
        diffs[:,j] = new_column
        diffs[j,:] = new_column
        diffs = np.delete(diffs, i, axis=0)
        diffs = np.delete(diffs, i, axis=1)
        clusters[j] = [clusters[i], clusters[j]]
        del clusters[i]
        return diffs, clusters

    def upgma(traits: Any, objects: List[str], output_file: str=None) -> None:
        """
            traits: A trait matrix 
            objects: A list of objects on which UPGMA is to be performed
            (output_file): A file to which intermediate values are to be written
                If this argument isn't present, intermediate values aren't logged

            returns: Newick representation of final tree
        """

        diffs = UPGMA.diffMatrixBuilder(traits)
        print("Difference Matrix: \n", diffs, '\n')
        clusters = [[thing] for thing in objects]
        assert(diffs.shape[0] == diffs.shape[1])
        assert(len(clusters) == diffs.shape[0])

 
        if output_file is not None:
            UPGMA.outputMatrix(output_file, diffs, clusters)

        while(len(clusters) > 2):
            diffs, clusters = UPGMA.upgmaStep(diffs, clusters)
            UPGMA.outputMatrix(output_file, diffs, clusters)

        return UPGMA.newickFromListOfLists(clusters)


if __name__ == "__main__":
    """
        This is what's run when you run python upgma.py. It takes the traitmatrix
        from "traitmatrix.txt", assigns A, B, C, ... to each of the objects on which
        UPGMA is being done, and constructs a tree based on the UPGMA. It also writes the 
        difference matrices to console and to the file out.txt
    """
    #TODO let user use other trait files and other outputs (so take those files as input)
    traits = np.loadtxt("traitmatrix.txt")
    objects = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ"[:traits.shape[0]])
    open("out.txt",'w').close()

    newick = UPGMA.upgma(traits, objects, "out.txt")
    print("The constructed tree, in Newick form, is: \n", newick)
    tree = Tree(newick+";")
    tree.show()
    



        


        





    
