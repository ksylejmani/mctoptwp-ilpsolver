import sys
import math
from ortools.linear_solver import pywraplp

class POI:
    def __init__(self,x,y,T,S,O,C,b=0, attribute_constraint=list()) -> None:
        """Initialize a POI"""
        self.x=x
        self.y=y
        self.T=T
        self.S=S
        self.O=O
        self.C=C
        self.b=b
        self.attribute_constraint=attribute_constraint

class MCTOPTWP_DATA:
    def __init__(self,M,N,B,max_attribute_constraint,pattern_sequence,pois:list) -> None:
        self.M=M
        self.N=N
        self.B=B
        self.max_attribute_constraint=max_attribute_constraint
        self.pattern_sequence=pattern_sequence
        self.Tmax=pois[0].C
        self.K=len(self.max_attribute_constraint)
        self.pois=pois
        self.P=2**32 # A large constant

class MCTOPTWP_ILP:
    def __init__(self, instance_type, instance_name) -> None:
        """Initiate object"""
        self.mctoptwp=self.parse_instance(instance_type,instance_name)
        self.travel_time=self.calculate_travel_time()
        pass

    def parse_instance(self, instance_type, instance_name):
        """Read data from file"""  
        with open('instances/MCTOPP-'+instance_type+'/'+instance_name,'r')  as file_read:
            # Read number of tours, number of locations, budget
            text_line=file_read.readline()
            text_list=text_line.split(' ')
            M=int(text_list[0])
            N=int(text_list[1])
            B=int(text_list[2])

            # Read attribute constraints
            text_line=file_read.readline()
            max_attribute_constraint=[int(i) for i in text_line.split(' ')]

            # Read pattern sequence
            text_line=file_read.readline()
            pattern_sequence=list()
            for i in range(M):
                text_line=file_read.readline()
                current_pattern=[int(i) for i in text_line.split(' ')]
                pattern_sequence.append(current_pattern)

            pois=list()

            # Read starting point
            text_line=file_read.readline()
            text_list=text_line.split(' ')
            x=float(text_list[1])
            y=float(text_list[2])
            T=float(text_list[3])
            S=float(text_list[4])
            O=int(text_list[5])
            C=int(text_list[6])
            start_poi=POI(x,y,T,S,O,C)
            pois.append(start_poi)

            # Read data for all points
            for i in range(N):
                text_line=file_read.readline()
                text_list=text_line.split(' ')
                x=float(text_list[1])
                y=float(text_list[2])
                T=float(text_list[3])
                S=float(text_list[4])
                O=int(text_list[5])
                C=int(text_list[6])
                b=int(text_list[7])
                attribute_constraint=list()
                for i in range(len(max_attribute_constraint)):
                    attribute_constraint.append(int(text_list[8+i]))
                poi=POI(x,y,T,S,O,C,b,attribute_constraint)
                pois.append(poi)
            
            pois.append(start_poi)
            mctoptwp_data=MCTOPTWP_DATA(M,N,B,max_attribute_constraint,pattern_sequence,pois)
            return mctoptwp_data

    def calculate_travel_time(self):
        """Calculate trave time"""
        travel_time={}
        for i in range(len(self.mctoptwp.pois)):
            for j in range(i+1,len(self.mctoptwp.pois)):
                travel_time[(i,j)]=math.sqrt(math.pow(self.mctoptwp.pois[i].x-self.mctoptwp.pois[j].x,2)+
                                             math.pow(self.mctoptwp.pois[i].y-self.mctoptwp.pois[j].y,2))
        return travel_time

    def transform_tuple(self, t:tuple)->tuple:
        """Transform tuple indeces"""
        min_index=min(t)
        max_index=max(t)
        return (min_index,max_index)
    
    def create_model(self, computation_time):
        """"Create Google OR-Tools model"""
        solver=pywraplp.Solver.CreateSolver('BOP')

        # Decision variables
        x = {}
        for i in range(self.mctoptwp.N):
            for j in range(self.mctoptwp.N):
                if i != j:
                    for d in range(self.mctoptwp.M):
                        x[(i, j,d)] = solver.BoolVar(name='x' + str(i)+ ',' + str(j) + ',' + str(d))
        y={}
        for i in range(self.mctoptwp.N):
            for d in range(self.mctoptwp.M):
                y[(i,d)]= solver.BoolVar(name='y' + str(i) + ',' + str(d))
        
        s={}
        for i in range(self.mctoptwp.N):
            for d in range(self.mctoptwp.M):
                s[(i,d)]= solver.IntVar(lb=0,ub=2**32,name='s' + str(i) + ',' + str(d))        
        
        # Constraints
        solver.Add(sum(sum(x[(0, j,d)] 
                           for j in range(1,self.mctoptwp.N)) 
                           for d in range(self.mctoptwp.M)) == 
                    sum(sum(x[(i, self.mctoptwp.N-1,d)] 
                           for i in range(self.mctoptwp.N-1)) 
                           for d in range(self.mctoptwp.M)), 
                           name='start from location 1 and end at location N')

        for o in range(1, self.mctoptwp.N-1):
            for d in range(self.mctoptwp.M):
                s1=sum(x[(i,o,d)] for i in range(self.mctoptwp.N-1) if i!=o)
                s2= sum(x[(o,j,d)] for j in range(1,self.mctoptwp.N) if j!=o)
                solver.Add( s1==s2==y[(o,d)], name='connectivity of each tour')
        
        for i in range(self.mctoptwp.N):
            for j in range(self.mctoptwp.N):
                if j!=i:
                    for d in range(self.mctoptwp.M):
                        solver.Add(s[(i,d)]+self.mctoptwp.pois[i].T+
                                   self.travel_time[self.transform_tuple((i,j))]-
                                   s[(j,d)]<=
                                   self.mctoptwp.P*(1-x[(i,j,d)]),
                                   name='time line of each tour')
        
        for i in range(1,self.mctoptwp.N-1):
            solver.Add(sum(y[(i,d)] for d in range(self.mctoptwp.M))<=1,
                       name='ensure that every location is visited at most once')
        
        for i in range(self.mctoptwp.N):
            for d in range(self.mctoptwp.M):
                solver.Add(self.mctoptwp.pois[i].O<=s[(i,d)])
        
        for i in range(self.mctoptwp.N):
            for d in range(self.mctoptwp.M):
                solver.Add(s[(i,d)]<=self.mctoptwp.pois[i].C)
        
        for d in range(self.mctoptwp.M):
            solver.Add(sum(self.mctoptwp.pois[i].T*y[(i,d)]
                           + sum(self.travel_time[self.transform_tuple((i,j))]*x[(i,j,d)] 
                                                 for j in range(1,self.mctoptwp.N) if i!=j)
                                                 for i in range(self.mctoptwp.N-1))<=
                                                 self.mctoptwp.Tmax,
                                                 name='limit the time budget')
        for k in range(self.mctoptwp.K):
            for d in range(self.mctoptwp.M):
                solver.Add(sum(self.mctoptwp.pois[i].attribute_constraint[k]*y[(i,d)] 
                               for i in range(1,self.mctoptwp.N-1))<=self.mctoptwp.max_attribute_constraint[k],
                               name='avoid the violation of the K attribute constraints')
        
        # Set objective function
        solver.Maximize(sum(sum(self.mctoptwp.pois[i].S*y[(i,d)] 
                            for i in range(1,self.mctoptwp.N-1) ) for d in range(self.mctoptwp.M)))
        
        # Configure solver
        solver.SetTimeLimit(computation_time*1000)
        solver.EnableOutput()
        
        # Call solver
        status=solver.Solve()

        # Print result
        if status == pywraplp.Solver.OPTIMAL:
            print('Solution:')
            print('Optimal value =', solver.Objective().Value())
        else:
            print('The problem does not have an optimal solution.')
            print('Objective value =', solver.Objective().Value())
        
        # Print visited POIs
        for i in range(self.mctoptwp.N):
            for d in range(self.mctoptwp.M):
                    if y[(i, d)].solution_value()==1:
                        print(y[(i, d)].name(), ' = ',  y[(i, d)].solution_value()) 
        
        # Print POI sequence
        for i in range(self.mctoptwp.N-1):
            for j in range(1,self.mctoptwp.N):
                if i!=j:
                    for d in range(self.mctoptwp.M):
                        if x[(i,j, d)].solution_value()==1:
                            print(x[(i,j, d)].name(), ' = ',  x[(i,j, d)].solution_value()) 

def help_function():
    """Define help function"""
    print("The solver should be callde usign the command 'python mctoptwp_ilp_solver.py instance_type computation_time instance_name.txt'")
    print("Example: python mctoptwp_ilp_solver.py Cordeau Scenario_input.txt 300")

if __name__=="__main__":
    """Create an object and use parameters"""
    instance_type="Cordeau"
    instance_name="Scenario_input.txt"
    computation_time=300
    mi=MCTOPTWP_ILP(instance_type, instance_name)
    mi.create_model(computation_time)

    # arguments=sys.argv
    # if len(arguments)!=4:
    #     help_function()
    #     exit()
    # else:
    #     instance_type = arguments[1]
    #     instance_name = arguments[2]
    #     computation_time=int(arguments[3])
    #     mi=MCTOPTWP_ILP(instance_type, instance_name)
    #     print(mi.mctoptwp.M)
    #     print(mi.mctoptwp.N)
    #     print(mi.mctoptwp.B)
    #     print(mi.mctoptwp.max_attribute_constraint)
    #     print(mi.mctoptwp.pois)