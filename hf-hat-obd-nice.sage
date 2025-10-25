## Copyright (C) 2025 Cagatay Kutluhan, Gordana Matic, Jeremy Van Horn-Morris, Andy Wand 
## Contact: kutluhan@buffalo.edu
## Based on hf-hat
## Copyright (C) 2015 Sucharit Sarkar
## Source: https://github.com/sucharit/hf-hat
## Contact: sucharit@math.princeton.edu

## hf-hat-obd-nice is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of
## the License, or (at your option) any later version.

## hf-hat-obd-nice is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with hf-hat; see COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

from datetime import datetime
import itertools
sage.parallel.decorate.parallel(p_iter='fork', ncpus=None)

class NiceHeegaardDiagram():

    def __init__(self,boundary_intersections,num_pointed_regions,name=None):
        print(datetime.now())
        print('Initializing diagram...')
        #Basic initialization. (Will convert boundary_intersections to its lexicographically smallest form.)   
        self.name=str(name)
        self.boundary_intersections=[]
        for intersections_list in boundary_intersections:
            smallest_indices=[2*ind_i for ind_i,i in enumerate(intersections_list[0::2]) if i==min(intersections_list[0::2])]
            next_intersections=[intersections_list[ind_i+1] for ind_i in smallest_indices]
            smallest_index=next(ind_i for ind_i in smallest_indices if intersections_list[ind_i+1]==min(next_intersections))
            self.boundary_intersections.append(intersections_list[smallest_index:]+intersections_list[:smallest_index])

        self.regions=range(len(self.boundary_intersections))#the regions
        self.regions_un=range(len(self.boundary_intersections)-num_pointed_regions)#the unpointed regions
        self.intersections=list(range(1+max(flatten(self.boundary_intersections))))#the intersection points

        #Euler measures of the regions, times 2 (unpointed data stored as a matrix).
        self.euler_measures_2=[2-len(self.boundary_intersections[R])/2 for R in self.regions]
        self.euler_measures_2_un=vector(ZZ,len(self.regions_un),self.euler_measures_2[:len(self.regions_un)])
        if min(self.euler_measures_2_un)>-1:
            self.is_nice=True#is a nice diagram
        else:
            self.is_nice=False
            raise Exception("The Heegaard diagram is not nice.")
        #boundary_mat[i][j] is the coefficient of the boundary of (boundary of the i-th region restricted to alpha circles) at the j-th intersection point  (unpointed data stored as a matrix).
        self.boundary_mat=[[-self.boundary_intersections[R][0::2].count(p)+self.boundary_intersections[R][1::2].count(p) for p in self.intersections] for R in self.regions]
        self.boundary_mat_un=matrix(ZZ,len(self.regions_un),len(self.intersections),self.boundary_mat[:len(self.regions_un)])
        
        #point_measures_4[i][j] is the point measure of i-th region at the j-th point, times 4  (unpointed data stored as a matrix).
        self.point_measures_4=[[self.boundary_intersections[R].count(p) for p in self.intersections] for R in self.regions]
        self.point_measures_4_un=matrix(ZZ,len(self.regions_un),len(self.intersections),self.point_measures_4[:len(self.regions_un)])
        
        self.periodic_domains=kernel(self.boundary_mat_un)#the free Abelian group of periodic domains
        self.pd_basis=basis(self.periodic_domains)#basis in echelon form for periodic domains
        self.b1=len(self.pd_basis)
        self.chern_num=[]#Chern numbers for the canonical spinc structure
        for D in self.pd_basis:
            self.chern_num.append(D.dot_product(self.euler_measures_2_un)/2)
        if self.chern_num==[0]*self.b1:
            self.torsion_spinc=True
        else:
            self.torsion_spinc=False
            self.div=gcd(self.chern_num)
        
        #intersections_on_alphas[i] is the ordered list of intersections on alpha_i (according to the orientation of alpha_i). Similarly, for beta.
        self.intersections_on_alphas=[]
        while len(flatten(self.intersections_on_alphas))<len(self.intersections):
            start_p=next(p for p in self.intersections if p not in flatten(self.intersections_on_alphas))
            new_circle=[]
            curr_p=start_p
            while curr_p!=start_p or new_circle==[]:
                new_circle.append(curr_p)
                found_next_point=False
                regions_with_curr_point=[(R,ind_p) for R in self.regions for ind_p in range(len(self.boundary_intersections[R])) if curr_p == self.boundary_intersections[R][ind_p]]
                for (R,ind_p) in regions_with_curr_point:
                    if not found_next_point:
                        if ind_p%2==0 and self.boundary_intersections[R][ind_p+1] not in new_circle:
                            found_next_point=True
                            curr_p=self.boundary_intersections[R][ind_p+1]
                        elif ind_p%2==1 and self.boundary_intersections[R][ind_p-1] not in new_circle:
                            found_next_point=True
                            curr_p=self.boundary_intersections[R][ind_p-1]
                if not found_next_point:#must have completed the cycle
                    curr_p=start_p
            self.intersections_on_alphas.append(new_circle)
            
        self.alphas=range(len(self.intersections_on_alphas))#the alpha circles
        self.contact_intersections=[]
        for a in self.alphas:
            self.contact_intersections.append(self.intersections_on_alphas[a][0])
            
        self.intersections_on_betas=[]
        while len(flatten(self.intersections_on_betas))<len(self.intersections):
            start_p=next(p for p in self.contact_intersections if p not in flatten(self.intersections_on_betas))
            new_circle=[]
            curr_p=start_p
            while curr_p!=start_p or new_circle==[]:
                new_circle.append(curr_p)
                found_next_point=False
                regions_with_curr_point=[(R,ind_p) for R in self.regions for ind_p in range(len(self.boundary_intersections[R])) if curr_p == self.boundary_intersections[R][ind_p]]
                for (R,ind_p) in regions_with_curr_point:
                    if not found_next_point:
                        if ind_p%2==0 and self.boundary_intersections[R][ind_p-1] not in new_circle:
                            found_next_point=True
                            curr_p=self.boundary_intersections[R][ind_p-1]
                        elif ind_p%2==1 and self.boundary_intersections[R][ind_p+1-len(self.boundary_intersections[R])] not in new_circle:
                            found_next_point=True
                            curr_p=self.boundary_intersections[R][ind_p+1-len(self.boundary_intersections[R])]
                if not found_next_point:#must have completed the cycle
                    curr_p=start_p
            self.intersections_on_betas.append(new_circle)

        self.betas=range(len(self.intersections_on_betas))#the beta circles
        if sorted(flatten(self.intersections_on_alphas))!=self.intersections:
            raise Exception("Alpha circles don't contain all the intersections")
        if sorted(flatten(self.intersections_on_betas))!=self.intersections:
            raise Exception("Beta circles don't contain all the intersections")

        self.adjacent_region_list=[]#list of unpointed regions adjacent to contact intersection points
        for p in self.contact_intersections:
            for region in self.regions_un:
                if p in self.boundary_intersections[region] and region not in self.adjacent_region_list:
                    self.adjacent_region_list.append(region)
        
        N=len(self.contact_intersections)
        self.regions_by_points=[]#unpointed region lists by contact intersection points
        for i in range(len(self.alphas)):
            regions_by_points_i=[]
            for region in self.adjacent_region_list:
                if self.contact_intersections[i] in self.boundary_intersections[region]:
                    regions_by_points_i.append(region)
            self.regions_by_points.append(regions_by_points_i)
    
        self.region_pairs=[]
        for i in range(len(self.alphas)):
            if len(self.regions_by_points[i])==2:
                self.region_pairs.append(self.regions_by_points[i])
        
        self.indicator_list=[]#the list of periodic domains in pd_basis and their multiplicities in the regions from adjacent_region_list
        for pdom in self.pd_basis:
            for region in self.adjacent_region_list:
                self.indicator_list.append(((self.pd_basis).index(pdom),region,pdom[region]))

        m=len(self.adjacent_region_list)
        n=len(self.pd_basis)
        self.indicator_matrix=matrix(ZZ,m,n)#matrix whose columns are multiplicities of periodic domains in adjacent unpointed regions at all contact intersection points
        for i in range(m):
            for j in range(n):
                self.indicator_matrix[i,j]=self.pd_basis[j][self.adjacent_region_list[i]]

        self.intersection_incidence=[[next(a for a in self.alphas if p in self.intersections_on_alphas[a]),next(b for b in self.betas if p in self.intersections_on_betas[b])] for p in self.intersections]#the i-th intersection point lies in alpha_a and beta_b, and has alpha.beta intersection number n, where [a,b,n]=intersection_incidence[i]
        for p in self.intersections:
            try:
                [a,b]=self.intersection_incidence[p]

                if len(self.intersections_on_alphas[a])>2 and len(self.intersections_on_betas[b])>2:#if any of a or b has 2 or fewer intersections, their orientations are not yet well-defined.
                    [a_pind,b_pind]=[(self.intersections_on_alphas[a]).index(p),(self.intersections_on_betas[b]).index(p)]
                    (R,ind_p)=next((R,ind_p) for R in self.regions for ind_p in range(len(self.boundary_intersections[R])/2) if self.boundary_intersections[R][2*ind_p]==p)
                    prev_b=self.boundary_intersections[R][2*ind_p-1]
                    next_a=self.boundary_intersections[R][2*ind_p+1]
                    intersection_number=-1#if a were oriented from p to next_a, and b oriented from prev_b to p, then the intersection.
                    if self.intersections_on_alphas[a][a_pind-1]==next_a:#a is actually oriented oppositely
                        intersection_number=-intersection_number
                    if self.intersections_on_betas[b][b_pind-1]!=prev_b:#b is actually oriented oppositely
                        intersection_number=-intersection_number
                    (self.intersection_incidence[p]).append(intersection_number)
                    
                elif len(self.intersections_on_alphas[a])==1 or len(self.intersections_on_betas[b])==1:
                    (self.intersection_incidence[p]).append(1)
                else:
                    print("WARNING: The current implementation of intersection numbers doesn't work if some alpha or beta circle has exactly 2 intersections.")
                    raise Exception("Couldn't orient the alpha and beta circles")#Comment out exception depending on how we feel.
            except ValueError:#intersection number at p already determined
                pass
        self.intersection_matrix=matrix(ZZ,len(self.alphas),len(self.betas))#the (a,b) entry is the intersection number between alpha circle a and beta circle b
        for a in self.alphas:
            for b in self.betas:
                self.intersection_matrix[a,b]=sum([n for [i,j,n] in self.intersection_incidence if (a==i and b==j)])
        
        #Now some more error checking. (Definitely not a complete collection.) If diagram is wrong, some error could also be raised by other parts of the program.
        if vector(ZZ,len(self.regions),(len(self.regions))*[1,])*matrix(ZZ,len(self.regions),len(self.intersections),self.point_measures_4)!=vector(ZZ,len(self.intersections),(len(self.intersections))*[4,]):
            raise Exception("Every intersection should have four corners.")
        if len(self.alphas)!=len(self.betas):
            raise Exception("Require same number of alpha and beta circles.")

        print(datetime.now())
        print('Initializing chain complex variables...')

        #generators is the list of hf-generators, and generator_reps are their representatives. Each generator is represented as a tuple of intersections; the i-th point will be on alpha_i.
        self.generator_reps=[]
        permutations=Arrangements(self.alphas,len(self.alphas)).list()
        for permu in permutations:
            possibilities=[[p for p in self.intersections_on_alphas[a] if self.intersection_incidence[p][1]==permu[a]] for a in self.alphas]
            self.generator_reps+=list(cartesian_product_iterator(possibilities))
        self.generators=range(len(self.generator_reps))#the generator names

        #Some more variables that will not be initialized initially, since we need to solve matrix equations.
        self.SpinC="Not yet initialized"#SpinC[i] will be the spinc structure of the i-th generator (numbered arbitrarily starting at 0).
        self.SpinC_structures="Not yet initialized"#the spinc structures
        self.genSpinC="Not yet initialized"#genSpinC[i] is a list of generators indices living in the i-th spinc structure.
        self.gradings="Not yet initialized"#the relative/absolute Maslov gradings of the generators
        self.differentials=dict()#a dictionary of Maslov index 1 positive domains between pairs of generators
        self.complex_initialized=False#will be set to true once the chain complex is initialized
        self.canonical_spinc_initialized=False#will be set to true once generators in canonical spinc are sorted
        self.boundary_maps_initialized=False#will be set to true once boundary maps in contact spinc are computed
        #The chain complexes will also not be initialized initially.

        #Variables needed to solve for domains.
        k=(self.boundary_mat_un.transpose()).nrows() #the same as the number of intersection points
        l=(self.boundary_mat_un.transpose()).ncols() #the same as the number of unpointed regions
        P, L, U = (self.boundary_mat_un.transpose()).LU() #B=PLU with B=self.boundary_mat_un.transpose() NEED TO BE MADE FASTER
        A=P*L #kxk invertible matrix
        self.Ainv=A.inverse()
        self.zrlist=[i for i in range(k) if U.row(i) == zero_vector(ZZ,l)] #list of zero rows of U
        self.plist=U.pivots() #pivot columns of U
        self.prlist=U.pivot_rows() #pivot rows of U
        self.nonplist=[i for i in range(l) if i not in self.plist] #non-pivot columns of U
        self.nonprlist=[i for i in range(k) if i not in self.prlist] #non-pivot rows of U
        self.U_red=U.matrix_from_rows_and_columns(self.prlist,self.plist) #extract pivot rows and columns of U
        self.Ainv_test=(self.Ainv).matrix_from_rows(self.zrlist) #extract the rows of Ainv corresponding to the zero rows of U
        self.Ainv_red=(self.Ainv).delete_rows(self.nonprlist) #remove the rows of Ainv corresponding to the nonpivot rows of U
        self.U_red_inv=(self.U_red).inverse() #U_red is an invertible (k-z)x(k-z) matrix where z is the length of nonprlist
        self.solve_mat=self.U_red_inv*self.Ainv_red #(k-z)xk matrix
        self.red_em_2_un=matrix(self.euler_measures_2[:len(self.regions_un)])
        self.red_em_2_un=self.red_em_2_un.delete_columns(self.nonplist)
        self.red_em_2_un_vect=vector(self.red_em_2_un)
        self.red_pm_4_un=matrix(self.point_measures_4[:len(self.regions_un)])
        self.red_pm_4_un=self.red_pm_4_un.delete_rows(self.nonplist)        
        print(datetime.now())
        print('Initialization is complete.')
        
    def find_domain(self,initial,final):
        #Finds a domain D from the initial generator to the final generator; returns (Maslov_index(D),D). Called by sort_generators.
        target_vect=zero_vector(len(self.intersections))
        target_vect_abs=zero_vector(len(self.intersections))
        for p in list(self.generator_reps[final]):
            target_vect[p]+=1
            target_vect_abs[p]+=1
        for p in list(self.generator_reps[initial]):
            target_vect[p]-=1
            target_vect_abs[p]+=1
        if self.Ainv_test*target_vect==zero_vector(len(self.zrlist)):#checks if there is an answer over the rationals
            solution=self.solve_mat*target_vect
            if [i for i in range(len(solution)) if not (solution[i]).is_integer()]==[]:
                maslov_4=2*solution.dot_product(self.red_em_2_un_vect)+solution*self.red_pm_4_un*target_vect_abs
                if maslov_4%4!=0:
                    return None
                else:
                    maslov=maslov_4/4
                    solution_list=list(solution)
                    for i in range(len(self.nonplist)):
                        solution_list.insert(self.nonplist[i],0)
                    domain=vector(ZZ,len(self.regions_un),solution_list)
                    return (maslov,domain)
            else:
                return None
        else:
            return None
        
    def sort_generators(self):
        #Sorts generators into spinc and Maslov gradings while recording domains between pairs of generators.
        print(datetime.now())
        print('Sorting generators into SpinC and Maslov gradings...')
        if self.complex_initialized:
            return True
        self.QSpinC=[] #list whose ith entry is the list of Chern numbers for the spinc structure associated to the ith generator
        self.SpinC_structures=[] #list of spinc structures with non-empty sets of generators
        self.SpinC=dict() #self.SpinC[i] is the spinc structure associated to the ith generator as specified by its Chern numbers and its torsion component       
        self.domains_stored=dict() #self.domains_stored[i] is a pair (m,D) where m is the relative Maslov grading between the ith generator and the base generator in the same spinc structure
        self.SpinC[0]=(self.chern_num,0)
        self.domains_stored[0]=(0,zero_vector(len(self.regions_un)))
        self.base_gens=[]
        for gen in self.generators:
            chern_num=[] #Chern numbers for the spinc structure associated to gen
            gen_vect_abs=vector(ZZ,len(self.intersections))
            for p in self.intersections:
                gen_vect_abs[p]=(self.generator_reps[gen]).count(p)+(self.generator_reps[gen]).count(p)
            for D in self.pd_basis:
                chern_num.append(D.dot_product(self.euler_measures_2_un)/2+(D*self.point_measures_4_un*gen_vect_abs)/4)
            (self.QSpinC).append(chern_num)
            if [foo for foo in self.base_gens if self.find_domain(gen,foo)!=None]==[]:
                (self.base_gens).append(gen)
                if [foo for foo in self.base_gens if foo<gen and self.QSpinC[foo]==chern_num]==[]:
                    self.SpinC[gen]=(chern_num,0)
                else:
                    self.SpinC[gen]=(chern_num,max([self.SpinC[foo][1] for foo in self.base_gens if foo<gen and self.QSpinC[foo]==chern_num])+1)
                self.domains_stored[gen]=(0,zero_vector(len(self.regions_un)))
            else:
                base_gen=min([foo for foo in self.base_gens if self.find_domain(gen,foo)!=None])
                self.SpinC[gen]=self.SpinC[base_gen]
                self.domains_stored[gen]=self.find_domain(gen,base_gen)
            if self.SpinC[gen] not in self.SpinC_structures:
                (self.SpinC_structures).append(self.SpinC[gen])
            else:
                continue
        self.genSpinC=[[gen for gen in self.generators if self.SpinC[gen]==S] for S in self.SpinC_structures] #list whose ith entry is a list of generators lying in the ith spinc structure
        self.div=dict() #list of integers whose ith entry is the divisibility of the ith spinc structure
        for S in self.SpinC_structures:
            ind_S=(self.SpinC_structures).index(S)
            if S[0]==[]:
                self.div[ind_S]=0
            else:
                self.div[ind_S]=gcd(S[0])
        self.gradings=[]
        for S in self.SpinC_structures:
            ind_S=(self.SpinC_structures).index(S)
            base_gen=self.genSpinC[ind_S][0]
            if S[0]==[0]*self.b1: #if S is a torsion spinc structure
                abs_gr=dict()
                abs_gr[base_gen]=0
                for gen in self.genSpinC[ind_S]:
                    abs_gr[gen]=abs_gr[base_gen]+self.domains_stored[gen][0]
                (self.gradings).append(abs_gr)
            else:
                rel_gr=dict()
                rel_gr[base_gen]=0
                for gen in self.genSpinC[ind_S]:
                    rel_gr[gen]=(rel_gr[base_gen]+self.domains_stored[gen][0])%self.div[ind_S]
                (self.gradings).append(rel_gr)

        self.cx0=[gen for gen in self.genSpinC[0] if self.gradings[0][gen]==0]
        self.cx1=[gen for gen in self.genSpinC[0] if self.gradings[0][gen]==1]
        
        self.complex_initialized=True
        self.canonical_spinc_initialized=True
        print('Sorting is complete.')
        print(datetime.now())

    def find_diffs(self,initial,final):
        #Finds all differentials from the initial generator to the final generator; returns a list of (Maslov_index(D),D). 
        alphas_used=[]
        for i in self.alphas:
            if self.generator_reps[initial][i]!=self.generator_reps[final][i]:
                alphas_used.append(i)
        diffs=[]
        if len(alphas_used)<=2:
            if self.SpinC[initial]!=None and self.SpinC[initial]==self.SpinC[final]: 
                if self.b1==0:
                    (m1,D1)=self.domains_stored[initial]
                    (m2,D2)=self.domains_stored[final]
                    (m,D)=(m1-m2,D1-D2)
                    if m==1 and [i for i in self.regions_un if D[i]<0]==[]:
                        diffs.append((m,D))
                    else:
                        diffs=[]
                else:
                    (m1,D1)=self.domains_stored[initial]
                    (m2,D2)=self.domains_stored[final]
                    domain=D1-D2
                    reduced_domain=[]
                    for region in self.adjacent_region_list:
                        reduced_domain.append(domain[region])
                    target_abs=vector(ZZ,len(self.intersections))
                    for p in self.intersections:
                        target_abs[p]=(self.generator_reps[final]).count(p)+(self.generator_reps[initial]).count(p)
                    multiplicities=[]
                    multiplicities.append([0]*len(self.adjacent_region_list))
                    indicator_alphas=[]
                    for i in alphas_used:
                        if self.contact_intersections[i] in self.generator_reps[final] and self.contact_intersections[i] not in self.generator_reps[initial]:
                            indicator_alphas.append(i)
                    if indicator_alphas!=[]:
                        uniq_adjacent_region_list=list(self.adjacent_region_list)
                        flat_adjacent_region_list=flatten(self.region_pairs)
                        truncated_indicator_alphas=[]
                        for a in indicator_alphas:
                            region_a_0=self.region_pairs[a][0]
                            region_a_1=self.region_pairs[a][1]
                            if flat_adjacent_region_list.count(region_a_0)==1 and flat_adjacent_region_list.count(region_a_1)==1:
                                truncated_indicator_alphas.append(a)
                                uniq_adjacent_region_list.remove(region_a_0)
                                uniq_adjacent_region_list.remove(region_a_1)
                            elif flat_adjacent_region_list.count(region_a_0)!=1 and uniq_adjacent_region_list.count(region_a_0)==1:
                                truncated_indicator_alphas.append(a)
                                uniq_adjacent_region_list.remove(region_a_0)
                            elif flat_adjacent_region_list.count(region_a_1)!=1 and uniq_adjacent_region_list.count(region_a_1)==1:
                                truncated_indicator_alphas.append(a)
                                uniq_adjacent_region_list.remove(region_a_1)
                        mult_at_contact=[[1,0],[0,1]]
                        choices=sorted(cartesian_product([[0,1] for i in range(len(truncated_indicator_alphas))]))
                        for choice in choices:
                            mult=[0]*len(self.adjacent_region_list)
                            for a in indicator_alphas:
                                region_a_0=self.region_pairs[a][0]
                                region_a_1=self.region_pairs[a][1]
                                if flat_adjacent_region_list.count(region_a_0)==1 and flat_adjacent_region_list.count(region_a_1)==1:
                                    a_ind=truncated_indicator_alphas.index(a)
                                    mult[(self.adjacent_region_list).index(region_a_0)]=mult_at_contact[choice[a_ind]][0]
                                    mult[(self.adjacent_region_list).index(region_a_1)]=mult_at_contact[choice[a_ind]][1]
                                elif flat_adjacent_region_list.count(region_a_0)!=1:
                                    if a in truncated_indicator_alphas:
                                        a_ind=truncated_indicator_alphas.index(a)
                                        mult[(self.adjacent_region_list).index(region_a_0)]=mult_at_contact[choice[a_ind]][0]
                                    mult[(self.adjacent_region_list).index(region_a_1)]=1-mult[(self.adjacent_region_list).index(region_a_0)]
                                elif flat_adjacent_region_list.count(region_a_1)!=1:
                                    if a in truncated_indicator_alphas:
                                        a_ind=truncated_indicator_alphas.index(a)
                                        mult[(self.adjacent_region_list).index(region_a_1)]=mult_at_contact[choice[a_ind]][1]
                                    mult[(self.adjacent_region_list).index(region_a_0)]=1-mult[(self.adjacent_region_list).index(region_a_1)]
                                else:
                                    pass
                            if mult not in multiplicities:
                                multiplicities.append(mult)       
                    m=len(self.adjacent_region_list)
                    K=len(multiplicities)
                    mult_matrix=Matrix(ZZ,m,K)
                    for i in range(m):
                        for j in range(K):
                            mult_matrix[i,j]=multiplicities[j][i]
                    domain_matrix=Matrix(ZZ,len(self.adjacent_region_list),len(multiplicities))
                    for i in range(len(self.adjacent_region_list)):
                        for j in range(len(multiplicities)):
                            domain_matrix[i,j]=reduced_domain[i]
                    delta_matrix=mult_matrix-domain_matrix
                    for k in range(len(multiplicities)):
                        delta=delta_matrix.column(k)
                        try:    
                            coef=(self.indicator_matrix).solve_right(delta)
                        except:
                            pass
                        else:
                            pD=[0]*len(self.regions_un)
                            for j in range(len(self.regions_un)):
                                for i in range(len(self.pd_basis)):
                                    pD[j]+=coef[i]*self.pd_basis[i][j]
                            new_domain=[domain[j]+pD[j] for j in range(len(self.regions_un))]
                            new_domain=vector(ZZ,new_domain)
                            maslov_4=2*new_domain.dot_product(self.euler_measures_2_un)+new_domain.dot_product(self.point_measures_4_un*target_abs)
                            maslov=maslov_4/4
                            if maslov==1 and [a for a in new_domain if a<0]==[]:
                                diffs.append((maslov,new_domain))
            else:
                pass
        else:
            diffs=[]
        return diffs
 
    def contributes(self,initial,final):
        #Checks if there is a Maslov index 1 positive domain from the initial generator to the final generator.
        if initial!=final and len(self.find_diffs(initial,final))%2==1:
            return True
        else:
            return False
        
    def find_diffs_from(self,initial): 
        #Finds all generators to which the initial generator has a differential.
        target_gens=[]
        for i in [i for i in self.alphas if self.generator_reps[initial][i]!=self.contact_intersections[i]]:
            p0=self.generator_reps[initial][i]
            sigma_i=self.intersection_incidence[p0][1]
            for q0 in [q for q in self.intersections_on_alphas[i] if q!=p0 and self.intersection_incidence[q][1]==sigma_i]: #find bigons
                gen=list(self.generator_reps[initial])
                gen[i]=q0
                gen=tuple(gen)
                final=(self.generator_reps).index(gen)
                if self.contributes(initial,final):
                    target_gens.append(final)
                else:
                    pass
            list_i=[j for j in self.alphas if j>i and self.generator_reps[initial][j]!=self.contact_intersections[j]]
            for j in list_i: #find rectangles
                p0=self.generator_reps[initial][i]
                p1=self.generator_reps[initial][j]
                sigma_i=self.intersection_incidence[p0][1]
                sigma_j=self.intersection_incidence[p1][1]
                for q0 in [q for q in self.intersections_on_alphas[i] if q!=p0 and self.intersection_incidence[q][1]==sigma_j]:
                    for q1 in [q for q in self.intersections_on_alphas[j] if q!=p1 and self.intersection_incidence[q][1]==sigma_i]:
                        gen=list(self.generator_reps[initial])
                        gen[i]=q0
                        gen[j]=q1
                        gen=tuple(gen)
                        final=(self.generator_reps).index(gen)
                        if self.contributes(initial,final):
                            target_gens.append(final)
                        else:
                            pass
        return target_gens
    
    def find_diffs_to(self,final):
        #Finds all generators from which there is a differential to the final generator.
        source_gens=[]
        for i in [i for i in self.alphas]:
            q0=self.generator_reps[final][i]
            sigma_i=self.intersection_incidence[q0][1]
            for p0 in [p for p in self.intersections_on_alphas[i] if p!=q0 and self.intersection_incidence[p][1]==sigma_i]: #first find bigons
                gen=list(self.generator_reps[final])
                gen[i]=p0
                gen=tuple(gen)
                initial=(self.generator_reps).index(gen)
                if self.contributes(initial,final):
                    source_gens.append(initial)
                else:
                    pass
            list_i=[j for j in self.alphas if j>i]
            for j in list_i: #find rectangles
                q0=self.generator_reps[final][i]
                q1=self.generator_reps[final][j]
                sigma_i=self.intersection_incidence[q1][1]
                sigma_j=self.intersection_incidence[q0][1]
                for p0 in [p for p in self.intersections_on_alphas[i] if p!=q0 and self.intersection_incidence[p][1]==sigma_i]:
                    for p1 in [p for p in self.intersections_on_alphas[j] if p!=q1 and self.intersection_incidence[p][1]==sigma_j]:
                        gen=list(self.generator_reps[final])
                        gen[i]=p0
                        gen[j]=p1
                        gen=tuple(gen)
                        initial=(self.generator_reps).index(gen)
                        if self.contributes(initial,final):
                            source_gens.append(initial)
                        else:
                            pass
        return source_gens
    
    def print_differentials(self,S):
        #Returns and records all Heegaard Floer differentials in in the spinc structure S.
        if not self.complex_initialized:
            self.sort_generators()
        ft=open(str(self.name)+'_differentials_in_spinc_'+str(S)+'.txt','a')
        Spinc=self.SpinC_structures[S]
        if len(self.genSpinC[S])>1:
            ft.write("In SpinC structure "+str(repr(Spinc))+'\n')
            for initial in self.genSpinC[S]:
                for final in [gen for gen in self.genSpinC[S] if self.contributes(initial,gen)]:
                    ft.write(str(repr(initial))+" in (relative) grading "+str(repr(self.gradings[S][initial]))+" has non-zero differential to "+str(repr(final))+'\n')
        ft.close()
        return None
            
    def plot_complex(self,S,output_file=None):
        #Plots the chain complex in the spinc structure S.
        if not self.complex_initialized:
            self.sort_generators()
        diffs=[]
        colors={'1': "blue", '(?)': "red"}
        heights=self.gradings
        for initial in self.genSpinC[S]:
            for final in self.genSpinC[S]:
                if self.contributes(initial,final):
                    diffs.append((initial,final,'1'))
                else:
                    continue
        G = DiGraph(diffs, multiedges=True) 
        Gplot = plot(G, edge_labels=false, edge_colors=G._color_by_label(colors), layout='acyclic', edge_labels_background= 'transparent', figsize=[40, 20], vertex_size=100)
        if output_file==None:
            Gplot.show() 
        else:
            Gplot.save(output_file)
    
    def compute_homology(self,S):
        #Computes the hat version of Heegaard Floer homology in spinc structure S.
        if not self.complex_initialized:
            self.sort_generators()
        print(datetime.now())
        print('Computing the differential of the chain complex...')
        chain_complex=dict()#This will be chain complex in spinc structure S. chain_complex[gr] will be the matrix from grading gr to gr-1.
        gradings_in_S=sorted(list(Set([self.gradings[S][g] for g in self.genSpinC[S]])))#the gradings that can appear in this spinc structure.
        for gr in gradings_in_S:
            if (gr-1) in gradings_in_S:
                generators_higher=[g for g in self.genSpinC[S] if self.gradings[S][g]==gr]
                generators_lower=[g for g in self.genSpinC[S] if self.gradings[S][g]==gr-1]

                boundary_matrix=matrix(GF(2),len(generators_higher),len(generators_lower))
                for ind_g,g in enumerate(generators_higher):
                    for ind_h,h in enumerate(generators_lower):
                        if self.contributes(g,h):
                            boundary_matrix[ind_g,ind_h]=1
                chain_complex[gr]=boundary_matrix
        print(datetime.now())
        print('Computing homology...')
        homology=[]
        homology_in_S=dict()
        ranks_in_S=dict()
        for gr in chain_complex:
            ranks_in_S[gr]=(chain_complex[gr]).rank()
        gradings_in_S=sorted(list(Set([self.gradings[S][gen] for gen in self.genSpinC[S]])))
        for gr in gradings_in_S:
            homology_in_S[gr]=len([gen for gen in self.genSpinC[S] if self.gradings[S][gen]==gr])
            if (gr-1) in gradings_in_S:
                homology_in_S[gr]-=ranks_in_S[gr]
            if (gr+1) in gradings_in_S:
                homology_in_S[gr]-=ranks_in_S[gr+1]
        print(datetime.now())
        
        return homology_in_S
    
    def sort_canonical_spinc(self):
        #Find generators in Maslov gradings 0 and 1 relative to the contact generator in the canonical spinc structure.
        print(datetime.now())
        print('Searching for generators in Maslov gradings 0 and 1 relative to the contact generator in the canonical SpinC structure...')
        if self.complex_initialized:
            self.canonical_spinc_initialized=True
        else:
            self.SpinC=dict()
            self.domains_stored=dict()
            self.cx=[]
            self.cx0=[]
            self.cx1=[]
            self.SpinC_structures=[(self.chern_num,0)]
            self.genSpinC=[[]]
            for gen in self.generators:
                if self.find_domain(gen,0)!=None:
                    (self.cx).append(gen)
                    (m,D)=self.find_domain(gen,0)
                    self.domains_stored[gen]=(m,D)
                    self.SpinC[gen]=(self.chern_num,0)
                    if (self.torsion_spinc and m==0) or (not self.torsion_spinc and m%self.div==0):
                        (self.cx0).append(gen)   
                    elif (self.torsion_spinc and m==1) or (not self.torsion_spinc and m%self.div==1):
                        (self.cx1).append(gen)
                    else:
                        continue
                    self.genSpinC[0].append(gen)
                else:
                    self.SpinC[gen]=None
            self.canonical_spinc_initialized=True
        print('Search is complete.')
        print(datetime.now())

    def cycle_reps(self,g):
        #Determines permutation types of generators in the Heegaard diagram of an open book decomposition.
        gen=self.generator_reps[g]
        cycle_gen=[]
        while len(flatten(cycle_gen))<len(gen):
            start_p=next(p for p in gen if p not in flatten(cycle_gen))
            new_cycle=[]
            curr_p=start_p
            while curr_p!=start_p or new_cycle==[]:
                new_cycle.append(curr_p)
                found_next_point=False
                next_p=next(p for p in gen if self.intersection_incidence[curr_p][1]==self.intersection_incidence[p][0])
                found_next_point=True
                curr_p=next_p
            if not found_next_point:
                curr_p=start_p
            cycle_gen.append(new_cycle)
        return cycle_gen
    
    def cycle_diff(self,initial,final):
        #Computes the cycle difference between two generators.
        cycle_diff=len(self.cycle_reps(initial))-len(self.cycle_reps(final))
        return cycle_diff
    
    def jplus_gr(self,gen):
        #Computes J_plus grading of generators in the canonical spinc structure as long as self.torsion_spinc==True
        if self.torsion_spinc==True and gen in self.genSpinC[0]:
            (maslov,domain)=self.find_domain(gen,0)
            euler_measure_2=domain.dot_product(self.euler_measures_2_un)
            iota=maslov-euler_measure_2
            jplus=iota+self.cycle_diff(gen,0)
            return jplus/2
        else:
            print('No absolute J_+ grading is defined for this generator.')
    
    def compute_bmaps(self):
        #Computes the Heegaard Floer differential and the J_+=0 and J_+=2 differentials between Maslov gradings 1 and 0 relative to the contact generator in the canonical spinc structure.
        if not self.canonical_spinc_initialized and not self.complex_initialized:
            self.sort_canonical_spinc() 
        print(datetime.now())
        print('Computing boundary maps between Maslov gradings 1 and 0 relative to the contact generator in the canonical SpinC structure...')
        m=len(self.cx0)
        n=len(self.cx1)
        self.dhat=Matrix(GF(2),m,n)
        self.d_zero=Matrix(GF(2),m,n)
        self.d_one=Matrix(GF(2),m,n)
        for i in range(m):
            for j in range(n):
                if self.contributes(self.cx1[j],self.cx0[i]):
                    self.dhat[i,j]=1
                    if (self.cycle_diff(self.cx1[j],self.cx0[i])==0 or self.cycle_diff(self.cx1[j],self.cx0[i])==-1):
                        self.d_zero[i,j]=1
                    else:
                        self.d_zero[i,j]=0
                    if self.cycle_diff(self.cx1[j],self.cx0[i])==1:
                        self.d_one[i,j]=1
                    else:
                        self.d_one[i,j]=0
                else:
                    self.dhat[i,j]=0
        self.ker_dhat=(self.dhat).right_kernel()
        self.img_dhat=(self.dhat).column_space()
        self.ker_d_zero=(self.d_zero).right_kernel()
        self.img_d_zero=(self.d_zero).column_space()
        self.ker_d_one=(self.d_one).right_kernel()
        self.img_d_one=(self.d_one).column_space()
        self.boundary_maps_initialized=True
        print(datetime.now())
        return None
        
    def check_contact_class(self):
        #Determines if the EH class is zero
        if not self.canonical_spinc_initialized and not self.complex_initialized:
            self.sort_canonical_spinc()
        if not self.boundary_maps_initialized:
            self.compute_bmaps()
        print(datetime.now())
        A=self.dhat
        m=A.nrows()
        n=A.ncols()
        EH=zero_vector(GF(2),m)
        EH[0]=1
        if EH in self.img_dhat:
            X=A.solve_right(EH)
            chain=[self.cx1[i] for i in range(n) if X[i]!=0]
            print("EH class is zero via "+str(chain)+".")
        else:
            print("EH class is non-zero.")
        print(datetime.now())
    
    def compute_order(self):
        #Computes the spectral order for the Heegaard diagram
        if not self.canonical_spinc_initialized and not self.complex_initialized:
            self.sort_canonical_spinc()
        if not self.boundary_maps_initialized:
            self.compute_bmaps()
        print(datetime.now())
        EH=zero_vector(GF(2),len(self.cx0))
        EH[0]=1
        d_zero=self.d_zero
        d_one=self.d_one
        if EH in self.img_d_zero: #if EH is not in the image of d_zero
            print("Order is 0.")
        elif d_zero.row(0)==0: #if there are differentials to the contact generator 
            print("Order is infinite.") 
        else:
            #Next pass to quotients to make induced d_one injective
            K0=d_one.right_kernel() #kernel of d_one
            K0mat=matrix(GF(2),K0.basis()).transpose() #matrix with columns given by a basis for K0
            d_zero_K0=(d_zero*K0mat).column_space() #d_zero(K0)
            d_zero_K0mat=matrix(GF(2),d_zero_K0.basis()) #matrix with rows given by a basis for d_zero(K0)
            Q=identity_matrix(len(self.cx0)) #matrix representing quotient by d_zero(K0)
            for i in range(len(self.cx0)-d_zero_K0.dimension(),len(self.cx0)):
                Q[i,i]=0
            comp_basis=[]
            for i in [i for i in range(len(self.cx0)) if i not in d_zero_K0mat.pivots()]:
                v=zero_vector(GF(2),len(self.cx0))
                v[i]=1
                comp_basis.append(v)
            Bchange=(matrix(GF(2),comp_basis).stack(matrix(GF(2),d_zero_K0.basis()))).transpose()
            Q=Bchange*Q*(Bchange).inverse()
            Qd_one=Q*d_one #d_one composed with quotient by d_zero(K0)
            K1=Qd_one.right_kernel() 
            while K1!=K0: #iterate the previous process
                K0=K1
                K0mat=matrix(GF(2),K0.basis()).transpose()
                d_zero_K0=(d_zero*K0mat).column_space()
                d_zero_K0mat=matrix(GF(2),d_zero_K0.basis())
                Q=identity_matrix(len(self.cx0)) #matrix representing quotient by d_zero(K0)
                for i in range(len(self.cx0)-d_zero_K0.dimension(),len(self.cx0)):
                    Q[i,i]=0
                comp_basis=[]
                for i in [i for i in range(len(self.cx0)) if i not in d_zero_K0mat.pivots()]:
                    v=zero_vector(GF(2),len(self.cx0))
                    v[i]=1
                    comp_basis.append(v)
                Bchange=(matrix(GF(2),comp_basis).stack(matrix(GF(2),d_zero_K0.basis()))).transpose()
                Q=Bchange*Q*(Bchange).inverse()
                Qd_one=Q*d_one
                K1=Qd_one.right_kernel()
            K=K1 #minimal subspace K such that the map induced by d_one from C1/K to C0/d_zero(K) is injective
            if K.dimension()==len(self.cx1):
                print("Order is infinite.")
            elif K==(Q*d_zero).right_kernel():
                print("Order is infinite.")
            else:
                Kmat=matrix(GF(2),K.basis())
                d_zero_K=(d_zero*Kmat.transpose()).column_space()
                d_zero_Kmat=matrix(GF(2),d_zero_K.basis())
                P=identity_matrix(len(self.cx1)).delete_columns([i for i in range(len(self.cx1)-K.dimension(),len(self.cx1))])
                comp_basis1=[]
                for i in [i for i in range(len(self.cx1)) if i not in Kmat.pivots()]:
                    v=zero_vector(GF(2),len(self.cx1))
                    v[i]=1
                    comp_basis1.append(v)
                Bchange1=(matrix(GF(2),comp_basis1).stack(Kmat)).transpose()
                Pred=Bchange1*P
                Q=identity_matrix(len(self.cx0)).delete_rows([i for i in range(len(self.cx0)-d_zero_K.dimension(),len(self.cx0))])
                comp_basis0=[]
                for i in [i for i in range(len(self.cx0)) if i not in d_zero_Kmat.pivots()]:
                    v=zero_vector(GF(2),len(self.cx0))
                    v[i]=1
                    comp_basis0.append(v)
                Bchange0=(matrix(GF(2),comp_basis0).stack(matrix(GF(2),d_zero_K.basis()))).transpose()
                Qred=Q*(Bchange0).inverse()
                d_zero_red=Qred*d_zero*Pred
                d_one_red=Qred*d_one*Pred

                EHmat=matrix(GF(2),Bchange0.inverse()*EH)
                EHmatred=EHmat.matrix_from_columns([j for j in range(0,len(self.cx0)-d_zero_K.dimension())])
                EHred=vector(GF(2),EHmatred.list())
                #Finally, compute the map delta and use it to compute order
                Id=identity_matrix(GF(2),d_one_red.ncols())
                for i in range(Id.nrows()):
                    row_i=d_one_red.solve_left(Id.row(i))
                    if i==0:
                        d_one_red_inv=matrix(GF(2),row_i)
                    else:
                        d_one_red_inv=d_one_red_inv.stack(row_i)
                delta=d_zero_red*d_one_red_inv
                V=d_zero_red.column_space()
                U0=d_one_red.column_space()
                delta_map=linear_transformation(delta)
                U1=(delta_map.inverse_image(U0)).intersection(U0)
                if U1==U0:
                    U=U1
                else:
                    while U1!=U0:
                        U0=U1
                        U1=(delta_map.inverse_image(U0)).intersection(U0)
                    U=U1
                W=(delta.right_kernel()).intersection(U)
                if EHred in span(V.basis()+W.basis(),GF(2)):
                    order=1
                else:
                    k=1
                    while EHred not in span(V.basis()+W.basis(),GF(2)):
                        k+=1
                        if W==((delta^k).right_kernel()).intersection(U):
                            order='infinite'
                            break
                        else:
                            W=((delta^k).right_kernel()).intersection(U)
                            if EHred in span(V.basis()+W.basis(),GF(2)):
                                order=k
                            else:
                                continue  
                print("Order is "+str(order)+".")
        print(datetime.now())