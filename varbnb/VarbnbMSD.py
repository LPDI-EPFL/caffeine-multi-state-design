import ipdb
import sys
from varbnb.SequenceNode import *
from varbnb.SequenceExpansionQueue import *
from ematrix.EMatrix import *
from dynamicAS.DynamicAS import *
from bp.BPnegativeStates import *
from mplp.MPLP import *
from config import config
class VarbnbMSD:
    ## allowedAAsPerMutRes: dictionary indexed by res_id_tups of the allowed AAs per mutable residue
    # state: the dictionary with the positive/negative states.
    def __init__(self, state, logfilename):
        self.logfile = open(logfilename, 'w')
        self.constraintMatrices = []
        self.state = state
        self.mergePositiveStateMatrices()
        self.constraintThresholdE = config.CONSTRAINT_THRESHOLD_E
        self.constraintThresholdDeltaE = config.CONSTRAINT_THRESHOLD_DELTA_E
    
    def log(self, s):
        if s is str:
            self.logfile.write(s)
        else:
            self.logfile.write(`s`)
        self.logfile.write("\n")
        self.logfile.flush()

    ## add a positive state to the problem. 
    def mergePositiveStateMatrices(self):
        # Merge positive state matrices. In the merged matrix, the rotamers of the same amino acid type 
        #     in residues that are common between the merged states are inserted as tuples.
        #    For example, if amino acid ASP at position ('A', 100) has 2 rotamers in each state, 
        #     then the merged matrix will have four states: 
        #    [('A', 100)][ASP][(0,0)], [('A', 100)][ASP][(0,1)], [('A', 100)][ASP][(1,0)], [('A', 100)][ASP][(1,1)].
        #    The energy of each "merged" rotamer is the sum of the previous one.
        self.log("Adding positive state; merging matrices ")
        assert(len(self.state['positive'].keys()) <= 2)

        if len(self.state['positive'].keys()) == 2:
            key1 = sorted(self.state['positive'].keys())[0] 
            key2 = sorted(self.state['positive'].keys())[1]
            self.positiveStatesMergedMatrix = EMatrix(self.state['positive'][key1]['ematrix'], self.state['positive'][key2]['ematrix'])

        else:
            key = sorted(self.state['positive'].keys())[0]
            self.positiveStatesMergedMatrix    = self.state['positive'][key]['ematrix']


        for statename in self.state['positive']:
            self.mergedPositiveAllowedAAs = {}
            for res_id_tup in self.positiveStatesMergedMatrix.getSortedResIds():
                self.mergedPositiveAllowedAAs[res_id_tup] = self.positiveStatesMergedMatrix.getAllowedAAs(res_id_tup)

    # This version assumes that two residues with the same number and chain (e.g. ('A', 100)) correspond to the same residue. 
    def optimize(self, initialSequence={}):

        # Compute constraint reference. 
        # TODO: This should be better implemented; it should not be tied to specific names.
        # Compute and add exact energy for reference states
        referenceE ={}
        referenceE['DeltaE'] = 0.0
        for statename in self.state['reference']:
            dynas = DynamicAS(self.state['reference'][statename]['ematrix'])
            stateTopNode = dynas.doAS()            
            myid = self.state['reference'][statename]['id']
            referenceE['DeltaE'] += config.SIGN[myid]*stateTopNode.fScore
            referenceE[myid] = {'name': statename, 'E': stateTopNode.fScore}
            self.log('Reference energy for state: {} id: {} = {:.2f}'.format(statename, myid, stateTopNode.fScore))
        self.log('Reference deltaE: {:.3f}'.format(referenceE['DeltaE']))
            


        # Start the A* SEQUENCE search. 
        queue = SequenceExpansionQueue()

        # Initialize empty root node.
        rootNode = SequenceNode(initialSequence)

        # Compute the sequence-expanded BP for the negative states, given the partial sequence of the node.
        # This step will compute a rotamer at every amino acid position that is considered the "best" rotamer for that position. 
        # This will be done for every negative state. 
        subtrahend = {}
        for statename in self.state['negative']:
            myemat = self.state['negative'][statename]['ematrix']
            bpNegSolver = BPnegativeStatesSolver(myemat, rootNode.partialSequence, statename)
            bpNegSolver.optimize()
            subtrahend[statename] = {}
            subtrahend[statename]['intra'], subtrahend[statename]['pair'] = bpNegSolver.computeSubtrahend()
            #self.log("Computing BP for negative state: {} energy: {}".format(statename, )
                    
        # Compute a bound on the positive/negative best score.
        mplpOptimizer = MPLP(self.positiveStatesMergedMatrix, subtrahend)
        energy, energyContributionPerResidueWithNegativeState = mplpOptimizer.optimizeEMPLP()
        rootNode.fScore = energy

        # Compute next residue to dynamically expand.
        # Compute fScore using MPLP again, but this time with no subtrahend matrices (i.e. only the positive state)
        mplpOptimizer = MPLP(self.positiveStatesMergedMatrix)
        energy, energyContributionPerResidueOnlyPositive = mplpOptimizer.optimizeEMPLP()
        # The next best residue is the one that has the highest differential between the two runs.
        rootNode.next_res_id_tup = self.computeNextBestRes(rootNode, energyContributionPerResidueOnlyPositive, energyContributionPerResidueWithNegativeState)

        # Place root node into queue.
        queue.insert(rootNode)

        # get the top element in the queue. 
        while (not queue.is_empty()):
            self.log("Entering loop")
            topNode = queue.pop()[1]
            self.log("Expanding node with fScore "+`topNode.fScore`+" and sequence ")
            self.log(topNode.partialSequence)

            # The first thing to check is that the node meets the constraints. 
            # Compute the constraint scores; if it is below the threshold, discard this node.
            discardNode = False
            # Compute constraint reference. 
            # Compute upper bounds on energy of synthetic bcl, synthetic binder and synthetic pair. 
            # TODO: The full constraints code should be reduced to a couple of lines.
            for statename in self.state['constraint']:
                myid = self.state['constraint'][statename]['id']
                mplpConstraintOptimizer = MPLP(self.state['constraint'][statename]['ematrix'])
                energy, dump = mplpConstraintOptimizer.optimizeEMPLP(topNode.partialSequence, subtrahend)
                if referenceE[myid]['E'] + self.constraintThresholdE < energy: 
                    self.log("Discarding node "+`topNode.partialSequence`+": energy bound: "+`energy`+" for    "+`statename`+\
                            " greater than wildtype + threshold: "+`(referenceE[myid]['E']+self.constraintThresholdE)`)
                    discardNode = True

            # Compute bound on delta E using MPLP/BP 
            subtrahend = {}
            for statename in self.state['constraint']:
                myid = self.state['constraint'][statename]['id']
                if config.SIGN[myid] == -1:
                    bpNegSolver = BPnegativeStatesSolver(self.state['constraint'][statename]['ematrix'], topNode.partialSequence)
                    bpNegSolver.optimize()
                    subtrahend[statename] = {}
                    subtrahend[statename]['intra'], subtrahend[statename]['pair'] = bpNegSolver.computeSubtrahend()
                    

            deltaE = None
            for statename in self.state['constraint']:
                myid = self.state['constraint'][statename]['id']
                if config.SIGN[myid] == 1:
                    mplpConstraintOptimizer = MPLP(self.state['constraint'][statename]['ematrix'])
                    deltaE, dump = mplpConstraintOptimizer.optimizeEMPLP(topNode.partialSequence, subtrahend)
            
            if deltaE and 'DeltaE' in referenceE and referenceE['DeltaE'] + self.constraintThresholdDeltaE < deltaE:
                self.log("Discarding node "+`topNode.partialSequence`+": energy bound: "+`deltaE`+" for deltaE"+\
                            " greater than wildtype + threshold: "+`referenceE['DeltaE'] + self.constraintThresholdDeltaE`)
                discardNode = True

            # If node's score is exact and we pulled it out of the queue it means that all other nodes are worse. we can stop. 
            if topNode.isExactScore and not discardNode:
                self.log("Found optimal sequence with energy: "+`topNode.fScore`)
                multistate_results = {}
                multistate_results['energy'] = topNode.fScore
                for statetype in topNode.stateNodes:
                    self.log(statetype+" states:")
                    multistate_results[statetype] = {}
                    for statename in topNode.stateNodes[statetype]:
                        self.log("Energy of "+statename+": "+`topNode.stateNodes[statetype][statename].fScore`)
                        multistate_results[statetype][statename] = {}
                        multistate_results[statetype][statename]['pdbfile'] = self.state[statetype][statename]['pdbfile']
                        multistate_results[statetype][statename]['energy'] = topNode.stateNodes[statetype][statename].fScore
                        multistate_results[statetype][statename]['residues'] = {}
                        for res_id_tup in topNode.stateNodes[statetype][statename].myConf:
                            multistate_results[statetype][statename]['residues'][`res_id_tup`] = {}
                            aa = topNode.stateNodes[statetype][statename].myConf[res_id_tup].keys()[0]
                            rot = topNode.stateNodes[statetype][statename].myConf[res_id_tup][aa]
                            multistate_results[statetype][statename]['residues'][`res_id_tup`]['aa'] = aa
                            multistate_results[statetype][statename]['residues'][`res_id_tup`]['Dih'] = \
                                    self.state[statetype][statename]['ematrix'].getDihedrals(res_id_tup, aa, rot)

                self.log(topNode.partialSequence)
                self.log("DeltaE of new sequence: "+`topNode.fScore`)
                if 'Delta' in referenceE:
                    self.log(" deltaE of referenceE:"+`referenceE['DeltaE']`+"; DeltaE of constraint: "+`topNode.constraintDeltaE`)
                return multistate_results

            # Check if the node has a fully assigned sequence by comparing to the allowed keys in the positive state. 
            # TODO: Right now we assume that the positive state has all the mutations allowed.
            if len(topNode.partialSequence.keys()) == len(self.mergedPositiveAllowedAAs) and not discardNode:
                topNode.fScore = 0.0
                topNode.constraintDeltaE = 0.0
                topNode.stateNodes = {}
                c ={'positive' : 1.0, 'negative' : -1.0}
                # Hack -- discard if positive is better than negative. 
                positive_scores = 0 


                # Compute and add exact energy for positive and negative states.
                for statetype in ['positive', 'negative']:
                    topNode.stateNodes[statetype] = {}
                    for statename in self.state[statetype]:
                        if 'homo' in self.state[statetype][statename]:
                            myPartialSeq = topNode.getHomoDimericSequence(self.state[statetype][statename]['homo'][0], self.state[statetype][statename]['homo'][1])
                            dynas = DynamicAS(self.state[statetype][statename]['ematrix'], myPartialSeq)
                        else:
                            dynas = DynamicAS(self.state[statetype][statename]['ematrix'], topNode.partialSequence)
                        stateTopNode = dynas.doAS()
                        topNode.fScore += c[statetype]*stateTopNode.fScore
                        topNode.stateNodes[statetype][statename] = stateTopNode
                        if statetype == 'positive': 
                            positive_scores = stateTopNode.fScore
                        elif statetype == 'negative':
                            if positive_scores > stateTopNode.fScore:
                                discardNode = True
                                self.log("Discarding node "+`topNode.partialSequence`+": negative state better than positive.")

                for statename in self.state['constraint']:
                    dynas = DynamicAS(self.state['constraint'][statename]['ematrix'], topNode.partialSequence)
                    myid = self.state['constraint'][statename]['id']
                    stateTopNode = dynas.doAS()
                    topNode.constraintDeltaE += config.SIGN[myid]*stateTopNode.fScore

                if not discardNode:
                    topNode.isExactScore = True
                    ## Reinsert into the queue.
                    queue.insert(topNode)

            elif (not discardNode):
                self.log("Residues in sequence: ")
                self.log(topNode.partialSequence)

                # Expand the next level of the node, according to node.next_res_id_tup
                next_res_id_tup = topNode.next_res_id_tup

                # For each amino acid allowed at the next residue of the node we will expand a new node:
                for aa in self.mergedPositiveAllowedAAs[next_res_id_tup]:
                    myResIdTup = next_res_id_tup
                    newNode = SequenceNode(parentPartialSequence=topNode.partialSequence, new_res_id_tup=myResIdTup, newAA=aa)
                    self.log(newNode.partialSequence)
                    
                    # Compute the sequence-expanded BP for the negative states, given the partial sequence of the node.
                    # This step will compute a rotamer at every amino acid position that is considered the "best" rotamer for that position. 
                    subtrahend = {}
                    for statename in self.state['negative']:
                        if 'homo' in self.state['negative'][statename]:
                            myPartialSeq = newNode.getHomoDimericSequence(self.state['negative'][statename]['homo'][0], self.state['negative'][statename]['homo'][1])
                            bpNegSolver = BPnegativeStatesSolver(self.state['negative'][statename]['ematrix'], myPartialSeq)
                        else: 
                            bpNegSolver = BPnegativeStatesSolver(self.state['negative'][statename]['ematrix'], newNode.partialSequence)
                        bpNegSolver.optimize()
                        subtrahend[statename] = {}
                        subtrahend[statename]['intra'], subtrahend[statename]['pair'] = bpNegSolver.computeSubtrahend()

                    # Compute fScore using MPLP 
                    mplpOptimizer = MPLP(self.positiveStatesMergedMatrix, subtrahend)
                    energy, energyContributionPerResidueWithNegativeState = mplpOptimizer.optimizeEMPLP(newNode.partialSequence)
                    newNode.fScore = energy
                    self.log("fScore of new node: = "+`newNode.fScore`)

                    # Compute next residue to dynamically expand.
                    # Compute fScore using MPLP again, but this time with no subtrahend matrices (i.e. only the positive state)
                    mplpOptimizer = MPLP(self.positiveStatesMergedMatrix)
                    energy, energyContributionPerResidueOnlyPositive = mplpOptimizer.optimizeEMPLP(newNode.partialSequence)
                    self.log("fScore of only positive state: = "+`energy`)
                    
    
                    # tHE Next best residue is the one that has the highest differential between the two runs.
                    newNode.next_res_id_tup = self.computeNextBestRes(newNode, energyContributionPerResidueOnlyPositive, energyContributionPerResidueWithNegativeState)
    
                    assert(newNode.next_res_id_tup not in newNode.partialSequence)
                    # Push the expanded node into the queue.
                    queue.insert(newNode)


        self.log("VarBNB Error: no solution was found. Likely constraints were not met. Exiting....")
        sys.exit(1)

    
    def computeAvailableAAsForNode(self, node):
        availableAAs = {}
        for res_id_tup in self.mergedPositiveAllowedAAs.keys():
            if res_id_tup in node.partialSequence:
                availableAAs[res_id_tup] = node.partialSequence[res_id_tup]
            else:
                availableAAs[res_id_tup] = self.mergedPositiveAllowedAAs[res_id_tup]
        return availableAAs

    #Computes the next best residue to expand based on energy contributions per residue in an only positive and a positive with negative matrix.
    def computeNextBestRes(self, node, energyContributionPerResidueOnlyPositive, energyContributionPerResidueWithNegativeState):
        # The next best residue is the one that has the highest differential between the two runs.
        max_diff_res_id_tup = ""
        max_diff = float("-inf")
        for res_id_tup in energyContributionPerResidueOnlyPositive.keys():
            if res_id_tup not in node.partialSequence:
                diff = energyContributionPerResidueOnlyPositive[res_id_tup] - energyContributionPerResidueWithNegativeState[res_id_tup]
                if diff > max_diff:
                    max_diff_res_id_tup = res_id_tup
                if max_diff_res_id_tup == "":
                    max_diff_res_id_tup = res_id_tup
        return max_diff_res_id_tup 
                 
