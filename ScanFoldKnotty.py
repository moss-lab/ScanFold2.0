import subprocess
import time

# change if knotty binary is not installed in path
KNOTTY_LOCATION = "knotty"

class KnottyObject: 
    counter = 0
    '''
    variables:
        sequence: sequence (string) 
        structure: structure (string)
        native mfe: mfe (float)
    if structure or mfe aren't proviced, the class will run knotty to produce them
    '''
    def __init__(self, sequence:str, structure: str = None, mfe: float = None):
        self.sequence = sequence
        self.structure = structure
        self.mfe = mfe
        KnottyObject.counter += 1

        if self.structure == None or self.mfe == None:
            self.structure, self.mfe = self.runKnotty(sequence)

    def runKnotty(self, sequence: str):
        '''
        input: a sequence to run knotty on
        output: structure, mfe

        use this in place of RNA.fold_compound in get_frag_feature_list in ScanFoldFunctions.py
        similar to algo == "rnastructure", return centroid = "NA" and native_mfe = 0.0
        '''
        #TODO: this is temporary, change it in the release version 
        starttime = time.time()
        knotty_location = KNOTTY_LOCATION
        args = [knotty_location,sequence,'-ns']
        #knotty_output = subprocess.Popen(args,stdout=subprocess.PIPE)
        #line = knotty_output.stdout.readline()
        knotty_output = subprocess.run(args, check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        #lines = knotty_output.decode().split('\n')
        print(str(knotty_output.stdout)) 
        print(str(knotty_output.stdout.decode))
        lines = str(knotty_output.stdout.decode()).splitlines()
        print(lines)
        structure = lines[1].split()[1]
        mfe = float(lines[1].split()[2])
        print(f"mfe: {mfe}")
        print(f"structure: {structure}")
        print(f'current nucleotide: {KnottyObject.counter}')
        endtime = time.time()
        runtime = endtime - starttime
        return structure,mfe
