## Copyright (C) 2024 Cagatay Kutluhan, Gordana Matic, Jeremy Van Horn-Morris, Andy Wand
## Contact: kutluhan@buffalo.edu

## makenice is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of
## the License, or (at your option) any later version.

## makenice is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with hf-hat-obd; see COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

class HeegaardDiagram():
    
    def __init__(self,boundary_intersections,num_pointed_regions):
        #Basic initialization. (Will convert boundary_intersections to its lexicographically smallest form.)
        self.boundary_intersections=[]
        self.num_pointed_regions=num_pointed_regions
        for intersections_list in boundary_intersections:
            smallest_indices=[2*ind_i for ind_i,i in enumerate(intersections_list[0::2]) if i==min(intersections_list[0::2])]
            next_intersections=[intersections_list[ind_i+1] for ind_i in smallest_indices]
            smallest_index=next(ind_i for ind_i in smallest_indices if intersections_list[ind_i+1]==min(next_intersections))
            self.boundary_intersections.append(intersections_list[smallest_index:]+intersections_list[:smallest_index])

        self.regions=range(len(self.boundary_intersections))#the regions
        self.regions_un=range(len(self.boundary_intersections)-num_pointed_regions)#the unpointed regions
        self.intersections=list(range(1+max(flatten(self.boundary_intersections))))#the intersection points
        self.pointed_regions = []
        for region in self.boundary_intersections[len(self.boundary_intersections)-self.num_pointed_regions:]:
            (self.pointed_regions).append(region)
        self.distances = [self.pointed_regions] #list whose i-th item is distance i regions
        self.max_distance = 1
        self.nice_boundary_intersections = [] #new list of boundary intersections after the diagram is made nice

    def beta_arcs(self,regions):
        #returns a list of beta arcs associated to the boundary of the domain composed of the given list of regions
        beta_arcs = []
        for region in regions:
            n = int(len(region)/2)
            beta_arcs.append([region[len(region)-1], region[0]])
            for i in range(n-1):
                beta_arcs.append([region[2*i+1], region[2*i+2]])
        return beta_arcs
    
    def outer_beta_arcs(self,regions):
    #returns a list of beta arcs associated to the boundary of the domain composed of the given list of regions
        outer_beta_arcs = []
        for region in regions:
            n = int(len(region)/2)
            if [region[len(region)-1], region[0]] not in outer_beta_arcs and [region[0], region[-1]] not in outer_beta_arcs:
                outer_beta_arcs.append([region[len(region)-1], region[0]])
            for i in range(n-1):
                if [region[2*i+1], region[2*i+2]] not in outer_beta_arcs and [region[2*i+2], region[2*i+1]] not in outer_beta_arcs:
                    outer_beta_arcs.append([region[2*i+1], region[2*i+2]])
        return outer_beta_arcs

    def adjacency(self,region,test_regions):
        #returns an ordered list which contains the given region as its initial item followed by a beta arc along which it is adjacent to the test_regions
        test_arcs = self.outer_beta_arcs(test_regions)
        n = int(len(region)/2)
        adjacent = False
        if [region[0],region[len(region)-1]] in test_arcs:
            adjacent = True
            adjacency = (region,[len(region)-1,0])
        else:
            for i in range(n-1):
                if [region[2*i+2],region[2*i+1]] in test_arcs:
                    adjacent = True
                    reordered_region = region[2*i+2:] + region[:2*i+2] #reorders the points on the boundary of the region cyclically so that the distance is achieved along the beta arc connecting the last point to the first point on the list
                    adjacency = (reordered_region,[region[2*i+1],region[2*i+2]])
                else:
                    continue
        if adjacent == True:
            return adjacency
        else:
            return None

    def compute_distances(self):
        current_test_regions = self.boundary_intersections[len(self.boundary_intersections)-self.num_pointed_regions:]
        control = self.boundary_intersections[0:len(self.boundary_intersections)-self.num_pointed_regions]
        current_distance = 0
        while control != []:
            regions_in_distance_foo = []
            regions_to_be_tested = control[:]
            new_test_regions = current_test_regions[:]
            for region in regions_to_be_tested:
                if self.adjacency(region, current_test_regions) != None:
                    regions_in_distance_foo.append(self.adjacency(region,current_test_regions)[0])
                    new_test_regions.append(region)
                    control.remove(region)
            current_distance += 1
            (self.distances).append(regions_in_distance_foo)
            current_test_regions = new_test_regions[:]
        self.max_distance = current_distance

    def make_nice(self):
        N = max(self.intersections)
        max_distance = len(self.distances)-1
        for i in range(max_distance):
            d = max_distance-i
            dist_d_bad_regions = []
            total_badness = 0
            for region in self.distances[d]:
                if len(region) > 4:
                    badness = int(len(region)/2)-2
                    dist_d_bad_regions.append(region)
                    total_badness += badness
                else:
                    continue
            current_regions = dist_d_bad_regions[:]
            while total_badness != 0:
                min_badness = min([int(len(region)/2)-2 for region in current_regions if int(len(region)/2)-2 != 0])
                current_region = next(region for region in current_regions if int(len(region)/2)-2 == min_badness)
                beta_arc = [current_region[len(current_region)-1],current_region[0]]
                adjacent_region = next(region for region in self.distances[d-1] if (beta_arc == [region[0],region[len(region)-1]] or beta_arc in [[region[2*i+2],region[2*i+1]] for i in range(int(len(region)/2-1))]))
                adjacent_region_betas = []
                for i in range(len(adjacent_region)/2-1):
                    adjacent_region_betas.append([adjacent_region[2*i+1],adjacent_region[2*i+2]])
                adjacent_region_betas.append([adjacent_region[len(adjacent_region)-1],adjacent_region[0]])####
                alpha_arc = [current_region[2],current_region[3]]
                target_region = next(region for region in flatten(self.distances,max_level=1) if alpha_arc in [[region[2*i+1],region[2*i]] for i in range(int(len(region)/2))])
                target_region_alphas = []
                for i in range(len(target_region)/2):
                    target_region_alphas.append([target_region[2*i],target_region[2*i+1]])
                new_rect = [current_region[0],current_region[1],current_region[2],N+1] #new rectangle region resulting from breaking the current_region two regions
                (self.distances[d]).append(new_rect)
                current_region[0] = N+2 #modify the current_region while breaking the orginal into two regions
                del current_region[1:3] #modify the current_region while breaking the orginal into two regions
                d_target = next(i for i in range(self.max_distance+1) if target_region in self.distances[i])
                num_rect = 0
                if len(target_region) == 4 and d_target > d-1: #if the target_region is a rectangle with distance at least d
                    curr_alpha_arc = alpha_arc
                    curr_target_region = target_region[:]
                    d_curr_target = d_target
                    while len(curr_target_region) == 4 and d_curr_target > d-1: #while satisfied, continue finger pushing isotopy             
                        drect_op = next(i+1 for i in range(self.max_distance+1) if [curr_target_region[2],curr_target_region[1]] in self.beta_arcs(self.distances[i]))
                        curr_target_region_ind = (self.distances[d_curr_target]).index(curr_target_region)
                        if [curr_alpha_arc[1],curr_alpha_arc[0]] == [curr_target_region[0],curr_target_region[1]]: #if the distance measuring beta arc is to the left
                            curr_alpha_arc = [curr_target_region[2],curr_target_region[3]]
                            self.distances[d_curr_target][curr_target_region_ind] = [curr_target_region[0],N+2,N+4,curr_target_region[3]]#first of three new rectangles resulting from finger pushing isotopy through the target_region
                            if d_curr_target == self.max_distance:
                                (self.distances).append([])
                                self.max_distance += 1
                            (self.distances[d_curr_target+1]).append([N+2,N+1,N+3,N+4])#next one of three new rectangles resulting from finger pushing isotopy through the target_region
                            if d_curr_target+2 < drect_op:
                                if d_curr_target+2 > self.max_distance:
                                    (self.distances).append([])
                                    self.max_distance += 1
                                (self.distances[d_curr_target+2]).append([N+1,curr_target_region[1],curr_target_region[2],N+3])#last one of three new rectangles resulting from finger pushing isotopy through the target_region
                            else:
                                if drect_op > self.max_distance:
                                    (self.distances).append([])
                                    self.max_distance += 1
                                (self.distances[drect_op]).append([curr_target_region[2],N+3,N+1,curr_target_region[1]])
                        else: #if the distance measuring beta arc is to the right
                            curr_alpha_arc = [curr_target_region[0],curr_target_region[1]]
                            self.distances[d_curr_target][curr_target_region_ind] = [curr_target_region[0],N+3,N+1,curr_target_region[3]]#first of three new rectangles resulting from pushing finger through the target_region
                            if d_curr_target == self.max_distance:
                                (self.distances).append([])
                                self.max_distance += 1
                            (self.distances[d_curr_target+1]).append([N+3,N+4,N+2,N+1])#next one of three new rectangles resulting from pushing finger through the target_region
                            if d_curr_target+2 < drect_op:
                                if d_curr_target+2 > self.max_distance:
                                    (self.distances).append([])
                                    self.max_distance += 1
                                (self.distances[d_curr_target+2]).append([N+4,curr_target_region[1],curr_target_region[2],N+2])#last one of three new rectangles resulting from pushing finger through the target_region
                            else:
                                if drect_op > self.max_distance:
                                    (self.distances).append([])
                                    self.max_distance += 1
                                (self.distances[drect_op]).append([curr_target_region[2],N+2,N+4,curr_target_region[1]])
                        curr_target_region = next(region for region in flatten(self.distances,max_level=1) if curr_alpha_arc in [[region[2*i+1],region[2*i]] for i in range(int(len(region)/2))])
                        d_curr_target = next(i for i in range(self.max_distance+1) if curr_target_region in self.distances[i])
                        N += 2
                        num_rect += 1
                        continue
                else: #if the finger pushing isotopy has to stop
                    curr_alpha_arc = alpha_arc[:]
                    curr_target_region = target_region[:]
                    d_curr_target = d_target 
                curr_target_region_ind = (self.distances[d_curr_target]).index(curr_target_region)
                target_region_alphas = []
                for i in range(len(curr_target_region)/2):
                    target_region_alphas.append([curr_target_region[2*i],curr_target_region[2*i+1]])

                if adjacent_region == curr_target_region:
                    new_bigon = [N+2,N+1]
                    (self.distances[d]).append(new_bigon)
                    beta_ind = adjacent_region_betas.index([beta_arc[1],beta_arc[0]])
                    if beta_ind == len(adjacent_region)-1:####
                        adjacent_region + [N-2*num_rect+1,N-2*num_rect+2]####
                    else:####
                        adjacent_region[2*beta_ind+2:2*beta_ind+2] = [N-2*num_rect+1,N-2*num_rect+2] #modified adjacent region
                    adjacent_region_alphas = []
                    for i in range(len(adjacent_region)/2):
                        adjacent_region_alphas.append([adjacent_region[2*i],adjacent_region[2*i+1]])
                    alpha_ind = adjacent_region_alphas.index([curr_alpha_arc[1],curr_alpha_arc[0]])
                    adjacent_region[2*alpha_ind+1:2*alpha_ind+1] = [N+2,N+1]
                    total_badness -= 1
                    N += 2
                else: # if adjacent_region and curr_target_region are distinct
                    new_bigon = [N+2,N+1]
                    if d_curr_target == self.max_distance:
                        (self.distances).append([])
                        self.max_distance += 1
                    (self.distances[d_curr_target+1]).append(new_bigon)
                    beta_ind = adjacent_region_betas.index([beta_arc[1],beta_arc[0]])
                    if beta_ind == len(adjacent_region)-1:####
                        adjacent_region + [N-2*num_rect+1,N-2*num_rect+2]####
                    else:####
                        adjacent_region[2*beta_ind+2:2*beta_ind+2] = [N-2*num_rect+1,N-2*num_rect+2] #modified adjacent region
                    alpha_ind = target_region_alphas.index([curr_alpha_arc[1],curr_alpha_arc[0]])
                    if len(curr_target_region) == 2 or d_curr_target != d:
                        total_badness -= 1
                    self.distances[d_curr_target][curr_target_region_ind][2*alpha_ind+1:2*alpha_ind+1] = [N+2,N+1]
                    N += 2
                continue
            continue
        for i in range(1,len(self.distances)):
            for region in self.distances[i]:
                (self.nice_boundary_intersections).append(region)
        for region in self.distances[0]:
            (self.nice_boundary_intersections).append(region)