import sys

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
        self.pois=pois

class MCTOPTWP_ILP:
    def __init__(self, instance_type, instance_name) -> None:
        """Initiate object"""
        self.mctoptwp=self.parse_instance(instance_type,instance_name)
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

    def create_model(self, computation_time):
        """"Create Google OR-Tools"""
        pass

def help_function():
    """Define help function"""
    print("The solver should be callde usign the command 'python mctoptwp_ilp_solver.py instance_type computation_time instance_name.txt'")
    print("Example: python mctoptwp_ilp_solver.py Cordeau Scenario_input.txt 300")

if __name__=="__main__":
    """Create an object and use parameters"""
    # instance_type="Cordeau"
    # instance_type="Scenario_input.txt"
    # computation_time=300
    arguments=sys.argv
    if len(arguments)!=4:
        help_function()
        exit()
    else:
        instance_type = arguments[1]
        instance_name = arguments[2]
        computation_time=int(arguments[3])
        mi=MCTOPTWP_ILP(instance_type, instance_name)
        print(mi.mctoptwp.M)
        print(mi.mctoptwp.N)
        print(mi.mctoptwp.B)
        print(mi.mctoptwp.max_attribute_constraint)
        print(mi.mctoptwp.pois)

