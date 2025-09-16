## This file is part of spectral_order
## Copyright (C) 2022 Cagatay Kutluhan
## Copyright (C) 2022 Jeremy Van Horn-Morris 
## Copyright (C) 2022 Andy Wand 
## Contact: kutluhan@buffalo.edu
## Based on hf-hat
## Copyright (C) 2015 Sucharit Sarkar
## Source: https://github.com/sucharit/hf-hat
## Contact: sucharit@math.princeton.edu

## spectral_order is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of
## the License, or (at your option) any later version.

## spectral_order is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with spec_order; see COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

class heegaard_diagram():

    def __init__(self,boundary_intersections,num_pointed_regions,name=None):
    
        """
        Number all the regions (complement of alpha,beta), starting at
        0, so that all the unpointed regions come first. Number all
        the intersection points between alpha,beta, again starting at
        0 (the numbering is arbitrary).
        
        The boundary_intersections is a list whose i-th element is a list of
        intersection points that appear on the boundary of the i-th
        region, according to the boundary orientation, so that the
        part of the boundary joining the first two points lie on an
        alpha circle. (We are implicitly assuming that all regions are
        planar with one boundary (which can intersect itself).)

        The num_pointed_regions is the number of pointed regions
        (always the last few).
        """

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
        
        #point_measures_4[i][j] is the point measure of i-th region at the j-th point, times 4  (unpointed data stored as a matrix).
        self.point_measures_4=[[self.boundary_intersections[R].count(p) for p in self.intersections] for R in self.regions]
        self.point_measures_4_un=matrix(ZZ,len(self.regions_un),len(self.intersections),self.point_measures_4[:len(self.regions_un)])
        
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
        
    def duplicate_arc(self,k):
        """
        Inputs a Heegaard Diagram H and a number k indicating the index
        in H.alphas of the alpha curve we would like to duplicate and
        produces a new list of unpointed regions in the new Heegaard
        diagram with one extra basepoint. The new list of unpointed
        regions is returned as a list of boundary intersection points
        as in the definition of heegaard_diagram.
        """

        #First label the newly introduced intersetion points preserving the labeling of the original intersection points.
        contact_intersection=self.intersections_on_alphas[k][0]
        self_int_on_alpha=[]
        for p in self.intersections_on_alphas[k]:
            if p in self.intersections_on_betas[k]:
                self_int_on_alpha.append(p)
        ind_alpha=self_int_on_alpha.index(contact_intersection)
        self_int_on_alpha_k=self_int_on_alpha[ind_alpha:]+self_int_on_alpha[:ind_alpha]
        self_int_on_beta=[]
        for p in self.intersections_on_betas[k]:
            if p in self.intersections_on_alphas[k]:
                self_int_on_beta.append(p)
        ind_beta=self_int_on_beta.index(contact_intersection)
        self_int_on_beta_k=self_int_on_beta[ind_beta:]+self_int_on_beta[:ind_beta]
        num_alpha=len(self.intersections_on_alphas[k])
        num_beta=len(self.intersections_on_betas[k])
        num_self=len(self_int_on_alpha_k)
        num_int=len(self.intersections)
        new_int_on_alphas=[]
        new_int_on_betas=[]
        new_int_on_alphas_k=self.intersections_on_alphas[k][:]
        ind=(self.intersections_on_betas[k]).index(contact_intersection)
        new_int_on_betas_k=self.intersections_on_betas[k][ind:]+self.intersections_on_betas[k][:ind]
        for i in range(num_self):
            p=self_int_on_alpha_k[num_self-i-1]
            if self.intersection_incidence[p][2]==1:
                (new_int_on_alphas_k).insert((new_int_on_alphas_k).index(p),num_int+num_self-i-1)
            elif self.intersection_incidence[p][2]==-1:
                (new_int_on_alphas_k).insert((new_int_on_alphas_k).index(p)+1,num_int+num_self-i-1)
        new_int_on_alphas.append(new_int_on_alphas_k)
        new_alpha_int=range(num_int+num_self,num_int+num_alpha+(2*num_self))
        new_int_on_alphas.append(new_alpha_int)
        num_int+=num_alpha+(2*num_self)
        for i in range(num_self):
            p=self_int_on_beta_k[(num_self-i)%num_self]
            if i==0:
                (new_int_on_betas_k).append(num_int-1)
            elif i!=0 and self.intersection_incidence[p][2]==1:
                (new_int_on_betas_k).insert((self.intersections_on_betas[k]).index(p)+1,new_int_on_alphas[1][(new_int_on_alphas[0]).index(p)-1])
            elif i!=0 and self.intersection_incidence[p][2]==-1:
                (new_int_on_betas_k).insert((self.intersections_on_betas[k]).index(p),new_int_on_alphas[1][(new_int_on_alphas[0]).index(p)-1])
        (new_int_on_betas).append(new_int_on_betas_k)
        new_beta_int=[]
        for i in range(len(new_int_on_betas[0])):
            if i==0:
                new_beta_int.append(new_int_on_alphas[1][0])
            elif i==1:
                new_beta_int.append(new_int_on_alphas[0][1])
            elif i in range(2,len(new_int_on_betas[0])):
                p=new_int_on_betas[0][i-1]
                if p in new_int_on_alphas[0] and self.intersection_incidence[p][2]==1:
                    new_beta_int.append(new_int_on_alphas[0][(new_int_on_alphas[0]).index(p)-1])
                elif p in new_int_on_alphas[1] and self.intersection_incidence[new_int_on_alphas[0][(new_int_on_alphas[1]).index(p)+1]][2]==1:
                    new_beta_int.append(new_int_on_alphas[1][new_int_on_alphas[1].index(p)-1])
                elif p in new_int_on_alphas[0] and self.intersection_incidence[p][2]==-1:
                    new_beta_int.append(new_int_on_alphas[0][(new_int_on_alphas[0]).index(p)+1])
                elif p in new_int_on_alphas[1] and self.intersection_incidence[new_int_on_alphas[0][(new_int_on_alphas[1]).index(p)+1]][2]==-1:
                    new_beta_int.append(new_int_on_alphas[1][(new_int_on_alphas[1]).index(p)+1])
                else:
                    new_beta_int.append(num_int)
                    num_int+=1
        (new_int_on_betas).insert(1,new_beta_int)

        #Next, generate the new_boundary_intersections for all the regions in the new Heegaard diagram using the labeling of the intersection points from above.
        num_pointed_regions=len(self.regions)-len(self.regions_un)
        regions_list=self.boundary_intersections[:]
        new_regions_list=[]
        oa=len(self.intersections_on_alphas[k])
        ob=len(self.intersections_on_betas[k])
        na=len(new_int_on_alphas[0])
        nb=len(new_int_on_betas[0])
        clone_regions=[]
        for region in regions_list:
            clone_region=[]
            rl=len(region)
            for ind_p in range(rl):
                p=region[ind_p]
                if p in self_int_on_alpha_k:
                    inda_p=(self.intersections_on_alphas[k]).index(p)
                    new_inda_p=(new_int_on_alphas[0]).index(p)
                    indb_p=(self.intersections_on_betas[k]).index(p)
                    new_indb_p=(new_int_on_betas[0]).index(p)
                    if ind_p%2==0 and region[(ind_p+1)%rl]==self.intersections_on_alphas[k][(inda_p+1)%oa] and self.intersection_incidence[p][2]==-1:
                        clone_region.insert(ind_p,new_int_on_alphas[1][(new_inda_p)%na])
                    elif ind_p%2==0 and region[(ind_p+1)%rl]==self.intersections_on_alphas[k][(inda_p+1)%oa] and self.intersection_incidence[p][2]==1:
                        clone_region.insert(ind_p,new_int_on_alphas[1][(new_inda_p-1)%na])
                    elif ind_p%2==1 and region[(ind_p-1)%rl]==self.intersections_on_alphas[k][(inda_p-1)%oa] and self.intersection_incidence[p][2]==-1:
                        clone_region.insert(ind_p,new_int_on_alphas[1][(new_inda_p-1)%na])
                    elif ind_p%2==1 and region[(ind_p-1)%rl]==self.intersections_on_alphas[k][(inda_p-1)%oa] and self.intersection_incidence[p][2]==1:
                        clone_region.insert(ind_p,new_int_on_alphas[1][(new_inda_p-2)%na])
                    elif ind_p%2==0 and region[(ind_p+1)%rl]==self.intersections_on_alphas[k][(inda_p-1)%oa] and self.intersection_incidence[p][2]==-1:
                        clone_region.insert(ind_p,p)
                    elif ind_p%2==0 and region[(ind_p+1)%rl]==self.intersections_on_alphas[k][(inda_p-1)%oa] and self.intersection_incidence[p][2]==1:
                        clone_region.insert(ind_p,new_int_on_betas[1][(new_indb_p+1)%nb])
                    elif ind_p%2==1 and region[(ind_p-1)%rl]==self.intersections_on_alphas[k][(inda_p+1)%oa] and self.intersection_incidence[p][2]==-1:
                        clone_region.insert(ind_p,new_int_on_betas[1][(new_indb_p+1)%nb])
                    elif ind_p%2==1 and region[(ind_p-1)%rl]==self.intersections_on_alphas[k][(inda_p+1)%oa] and self.intersection_incidence[p][2]==1:
                        clone_region.insert(ind_p,p)
                elif p in self.intersections_on_alphas[k] and p not in self.intersections_on_betas[k]:
                    inda_p=(self.intersections_on_alphas[k]).index(p)
                    new_inda_p=(new_int_on_alphas[0]).index(p)
                    if ind_p%2==0 and region[(ind_p+1)%rl]==self.intersections_on_alphas[k][(inda_p+1)%oa]:
                        clone_region.insert(ind_p,new_int_on_alphas[1][(new_inda_p-1)%na])
                    elif ind_p%2==1 and region[(ind_p-1)%rl]==self.intersections_on_alphas[k][(inda_p-1)%oa]:
                        clone_region.insert(ind_p,new_int_on_alphas[1][(new_inda_p-1)%na])
                    elif ind_p%2==0 and region[(ind_p+1)%rl]==self.intersections_on_alphas[k][(inda_p-1)%oa]:
                        clone_region.insert(ind_p,p)
                    elif ind_p%2==1 and region[(ind_p-1)%rl]==self.intersections_on_alphas[k][(inda_p+1)%oa]:
                        clone_region.insert(ind_p,p)
                elif p not in self.intersections_on_alphas[k] and p in self.intersections_on_betas[k]:
                    indb_p=(self.intersections_on_betas[k]).index(p)
                    new_indb_p=(new_int_on_betas[0]).index(p)
                    if ind_p%2==0 and region[(ind_p-1)%rl]==self.intersections_on_betas[k][(indb_p-1)%ob]:
                        clone_region.insert(ind_p,new_int_on_betas[1][(new_indb_p+1)%nb])
                    elif ind_p%2==1 and region[(ind_p+1)%rl]==self.intersections_on_betas[k][(indb_p-1)%ob]:
                        clone_region.insert(ind_p,p)
                    elif ind_p%2==0 and region[(ind_p-1)%rl]==self.intersections_on_betas[k][(indb_p+1)%ob]:
                        clone_region.insert(ind_p,p)
                    elif ind_p%2==1 and region[(ind_p+1)%rl]==self.intersections_on_betas[k][(indb_p+1)%ob]:
                        clone_region.insert(ind_p,new_int_on_betas[1][(new_indb_p+1)%nb])
                else:
                    clone_region.insert(ind_p,p)
            clone_regions.append(clone_region)
        clone_regions_un=clone_regions[0:len(self.boundary_intersections)-num_pointed_regions]
        pointed_regions_list=clone_regions[len(self.boundary_intersections)-num_pointed_regions:]
        alpha_strip_regions=[]
        for j in range(na-1):
            alpha_strip_regions.append([new_int_on_alphas[0][j+1],new_int_on_alphas[0][(j+2)%na],new_int_on_alphas[1][j+1],new_int_on_alphas[1][j]])
        beta_strip_regions=[]
        for j in range(nb-1):
            test=[new_int_on_betas[1][j+1],new_int_on_betas[0][j],new_int_on_betas[0][j+1],new_int_on_betas[1][(j+2)%nb]]
            tl=len(test)
            test_list=[]
            for i in range(tl):
                new_test=test[2*i:]+test[:2*i]
                test_list.append(new_test)
            already=[]
            for i in range(len(test_list)):
                if test_list[i] in alpha_strip_regions+clone_regions_un:
                    already.append(test_list[i])
            if already==[]:
                beta_strip_regions.append([new_int_on_betas[1][j+1],new_int_on_betas[0][j],new_int_on_betas[0][j+1],new_int_on_betas[1][(j+2)%nb]])

        pointed_regions_list.append([new_int_on_alphas[0][0],new_int_on_alphas[0][1],new_int_on_alphas[1][0],new_int_on_alphas[1][na-1]])
        new_boundary_intersections=alpha_strip_regions[:1]+beta_strip_regions[:1]+alpha_strip_regions[1:]+beta_strip_regions[1:]+clone_regions_un+pointed_regions_list

        return new_boundary_intersections

def arc_duplicator(H,multiplicities):
    #Takes a Heegaard diagram, a list of sequential alphas starting from zero, and a list of corresponding desired multiplicities, and outputs a diagram (with the appropriate basepoints)
    HD=H
    arc_groups=[]
    for a in H.alphas:
        arcs_a=[a]
        if multiplicities[a]==1:
            HD=HD
        else:
            i=0
            while i<multiplicities[a]-1:
                num_basepoints=len(HD.regions)-len(HD.regions_un)
                if i==0:
                    boundary_intersections=HD.duplicate_arc(a)
                else:
                    boundary_intersections=HD.duplicate_arc(len(HD.alphas)-1)
                num_basepoints+=1
                HD=heegaard_diagram(boundary_intersections,num_basepoints,HD.name)
                arcs_a.append(len(HD.alphas)-1)
                i+=1
        arc_groups.append(arcs_a)
    suffix=""
    for a in H.alphas:
        suffix+=str(multiplicities[a])
    diagram=str(HD.name)+"_"+str(suffix)+'\n'+str(HD.boundary_intersections)+'\n'+'\n'
    intersections=str(HD.intersections_on_alphas)+'\n'+'\n'
    parallel_arcs=str(arc_groups)+'\n'+'\n'
    ft=open(str(HD.name)+'.txt','a')
    ft.write(diagram)
    ft.write(intersections)
    ft.write(parallel_arcs)
    ft.close
    return HD.boundary_intersections