proc cells_swapping {cellName_list Vt_type} {
    foreach cellName $cellName_list {
        set cell [get_cells $cellName]
        set ref_name [get_attribute $cell ref_name]
        set original_vt [get_attribute $cell lib_cell.threshold_voltage_group]

        if {$original_vt != $Vt_type} {
            set library_name "CORE65LP${Vt_type}"
            if {$original_vt == "LVT"} {
                regsub {_LL} $ref_name "_LH" new_ref_name
            } elseif {$original_vt == "HVT"} {
                regsub {_LH} $ref_name "_LL" new_ref_name
            }
            size_cell $cell "${library_name}/${new_ref_name}"
        }
    }
}

proc swapping_to_HVT {} {
    foreach_in_collection cell [get_cells] {
        set cell_name [get_attribute $cell full_name]
        set ref_name [get_attribute $cell ref_name]
        set library_name "CORE65LPHVT"
        # you can also extract the group from the name using regex
        regsub {_LL} $ref_name "_LH" new_ref_name
        size_cell $cell "${library_name}/${new_ref_name}"
    }
}

proc constraints_met {slackThreshold maxFanoutEndpointCost cell_list} {
    update_timing -full
    set ret_val 1
    # check slack
    set wrt_slack [get_attribute [get_timing_paths] slack]
    if {$wrt_slack < 0} {
        # Slack not met
        set ret_val 0
    }
    # check FEC
    foreach elem $cell_list {
        set list_cells_in_worst_path {}
        set cell [get_cells [lindex $elem 4]]
        # we consider only violating endopoints, which are the ones corresponding to inputs or outputs for which we have a path violating the slack constraint
        set paths [get_timing_paths -through $cell -nworst 1 -max_paths 10000 -slack_lesser_than $slackThreshold]
        set cell_fanout_endpoint_cost 0.0
        set worst_slack 100
        foreach_in_collection path $paths {
            set index 0
            set worst_found 0
            set curr_slack [get_attribute $path slack]
            # to determine if the path that we are considering is the worst
            if {$curr_slack < $worst_slack} {
                set worst_found 1
                set worst_slack $curr_slack
                set list_cells_in_worst_path {} 
            }
            # to determine the list of cells in the current worst path
            foreach_in_collection timing_point [get_attribute $path points] {
                # consider the cells only once, the list gives us both input and output pins for the same cell
                if {[expr {$index % 2}]==0} {
                    set point_object [get_attribute $timing_point object]
                    set point_name [get_attribute $point_object full_name]
                    regsub {/.*} $point_name "" point_name
                    set found [regexp {^N.*} $point_name]
                    if {$found == 0 && $worst_found == 1} {
                        if {[lsearch $list_cells_in_worst_path $point_name] == -1} {
                            # append a new cell to the list of cells in the current worst path
                            lappend list_cells_in_worst_path $point_name
                        }
                    }
                }
                set index [expr {$index + 1}]
            }
            set this_cost [expr $slackThreshold - [get_attribute $path slack]]
            set cell_fanout_endpoint_cost [expr $cell_fanout_endpoint_cost + $this_cost]
        }
        # ::current_FECs is a global variable
        lappend ::current_FECs $cell_fanout_endpoint_cost
        lappend ::cells_in_paths $list_cells_in_worst_path
        lappend ::worst_slacks $worst_slack
        # verify if the FEC is coherent with the constraints
        if {$cell_fanout_endpoint_cost >= $maxFanoutEndpointCost} {
            set cell_name [get_attribute $cell full_name]
            set cell_ref_name [get_attribute $cell ref_name]
            set ret_val 0
        }
    }
    return $ret_val
}
proc computeP {weight_deltaL deltaL FEC slack max_deltaL max_FEC min_slack iter} {

    # the metric below is the best (...that we found)
    if {$slack > 0} { 
        set ret_val [expr {($weight_deltaL**4*$max_deltaL/($deltaL))+(1.2*$FEC/$max_FEC)}] 
    } else {
        set ret_val [expr {-0.4*($slack/$min_slack)*($max_deltaL/($deltaL))**($weight_deltaL)+(1.2*$FEC/$max_FEC)*(1.1**$iter)}]

    }
    return $ret_val
}

proc dualVth {slackThreshold maxFanoutEndpointCost} {
    # nested list, structured as : {{ref_name, deltaL, FEC, p, full_name, cells_in_paths, worst_slack},{...}}
    set cell_list {}
    set index 0
    set leakage_initial [get_attribute [current_design] leakage_power]
    # compute some parameters
    foreach_in_collection cell [get_cells] {
        lappend cell_list {0 0 0 0 0 0 0}
        # for now set as deltaL only the LVT leakage
        lset cell_list $index 1 [get_attribute $cell leakage_power]
        # set the cell name as the element in the first place of the nested list
        lset cell_list $index 0 [get_attribute $cell ref_name]
        # set the cell name as the element in the first place of the nested list
        lset cell_list $index 4 [get_attribute $cell full_name]
        set index [expr {$index + 1}]
    }
    swapping_to_HVT
    set index 0
    # delta leakage computation
    foreach elem $cell_list {
        # extract the cell with the name from the list
        set cell [get_cells [lindex $elem 4]]
        set hvt_leakage [get_attribute $cell leakage_power]
        lset cell_list $index 1 [expr { [lindex $elem 1] - $hvt_leakage }]
        set index [expr {$index + 1}]
    }
    set index 0
    set ::current_FECs {}
    set ::cells_in_paths {}
    set ::worst_slacks {}
    # initial parameters computation inside the constraints_met call
    set initial_met [constraints_met $slackThreshold $maxFanoutEndpointCost $cell_list]
    # if constraints are met with all HVT cells return 1
    if {$initial_met == 1} {
        return 1
    }
    set index 0
    # assign the parameters computed to the corresponding elements of the cell_list nested list
    foreach cell $cell_list {
        lset cell_list $index 2 [lindex $::current_FECs $index]
        lset cell_list $index 5 [lindex $::cells_in_paths $index]
        lset cell_list $index 6 [lindex $::worst_slacks $index]
        set index [expr {$index+1}]  
    }
    # list for LVT cells in the current configuration
    set lvt_cells {}
    # list for HVT cells in the current configuration
    set hvt_cells {}
    # at least one working configuration found (in a single execution with fixed p params)
    set working_found 0
    set max_deltaL 0
    set max_FEC 0
    # Find max values for deltaL and FEC
    foreach elem $cell_list {
        set curr_deltaL [lindex $elem 1]
        if { $curr_deltaL > $max_deltaL} {
            set max_deltaL [lindex $elem 1]
        }
        if { [lindex $elem 2] > $max_FEC} {
            set max_FEC [lindex $elem 2]
        }
    }
    # counter of iterations of the binary search
    set iter 0
    # Binary search's step by which we have to move
    set step [expr {[llength $cell_list]/4}]
    # number of cells to change
    set to_change [expr {[llength $cell_list]/2}]
    # outcome of the previous constraints check
    set prev_outcome 0
    # at least one working configuration found
    set working_found 0
    # Used to evaluate the quality of the solutions
    set normalized_power 10
    # number of divisions of an odd step, useful to explore the whole array even after the binary search reaches step 0
    set odd_divisions 0
    # infinite loop, break under certain conditions (emulation of do while)
    while { 0 < 1 } {
        # It is not possible to proceed with the binary search
        if {$step == 0} {
            break
        }
        set index 0
        # Find max values again, but only on the HVT cells that are left
        set max_FEC 0
        set min_slack 0
        set max_ncells 0
        set max_deltaL 0
        foreach elem [lrange $cell_list [expr {($to_change - $step*2)-1}] [expr {[llength $cell_list]-1}]] {
            if { [lindex $elem 2] > $max_FEC} {
                set max_FEC [lindex $elem 2]
            }
            if { [lindex $elem 6] < $min_slack && [lindex $elem 6] < 0} {
                set min_slack [lindex $elem 6]
            }
            if { [lindex $elem 1] > $max_deltaL} {
                set max_deltaL [lindex $elem 1]
            }
        }
        set weight_deltaL 0
        set divisor 0
        set index 0
        # to compute the weight_deltaL to be used in the computation of p
        while {$index < [llength [lrange $cell_list [expr {($to_change - $step*2)-1}] [expr {[llength $cell_list]-1}]]]} {
            if {[lindex [lrange $cell_list [expr {($to_change - $step*2)-1}] [expr {[llength $cell_list]-1}]] $index 6] < 0} {
                # compute the numerator of the sum step by step
                set weight_deltaL [expr {$weight_deltaL - [lindex [lrange $cell_list [expr {($to_change - $step*2)-1}] [expr {[llength $cell_list]-1}]] $index 1] * [lindex [lrange $cell_list [expr {($to_change - $step*2)-1}] [expr {[llength $cell_list]-1}]] $index 6]}]
                # compute the denominator step by step
                set divisor [expr {$divisor - $max_deltaL * $min_slack}]
            }
            set index [expr {$index + 1}]
        }
        # to avoid divisions by zero
        if {$divisor == 0} {
            set divisor 1
        }
        set weight_deltaL [expr {$weight_deltaL/$divisor}]
        # To make the min slack positive
        set min_slack [expr {0 - $min_slack}]
        if {$max_FEC == 0} { ;# to avoid infinite p values
            set max_FEC 1
        }
        if {$min_slack == 0} { ;# to avoid infinite p values
            set min_slack 1
        }
        set index 0
        # Compute p for all the nodes
        foreach elem $cell_list {
            # update all the metrics used for the computation of p
            if {$prev_outcome == 2} {
                lset cell_list $index 2 [lindex $::current_FECs $index]
                lset cell_list $index 5 [lindex $::cells_in_paths $index]
                lset cell_list $index 6 [lindex $::worst_slacks $index]
            }
            lset cell_list $index 3 [computeP $weight_deltaL [lindex $cell_list $index 1] [lindex $cell_list $index 2] [lindex $cell_list $index 6] $max_deltaL $max_FEC $min_slack $iter]
            set index [expr {$index + 1}]
        }
        if {$prev_outcome == 2} {
            # Sort if constraints have not been met in the previous iteration
            # but not if you are at the first iteration
            if { [expr {$to_change - $step*2}] == 0} {
                set index 0
                # the following loops are used to penalize cells that are in the same path of cells selected for swapping to LVT
                # for each cell among the cells that are left in the list (the HVT ones)
                while {$index < $to_change} {
                    set curr_cell [lindex $cell_list $index 4]
                    set explore_path 0
                    set path_list [lindex $cell_list $index 5]
                    # for each cell in the worst path of the selected cell
                    while {$explore_path < [llength [lindex $cell_list $index 5]]} {
                        set scan_cell_list $index
                        # scan the cell_list array to find the selected cell among the ones in the most critical path passing through the cell selected at the beginning 
                        while {$scan_cell_list < [llength $cell_list]} {
                            if {[lindex $cell_list $scan_cell_list 4] == [lindex $path_list $explore_path]} {  
                                if {[lindex $cell_list $scan_cell_list 6] >= 0} {
                                    # penalize it a lot
                                    lset cell_list $scan_cell_list 3 [expr {[lindex $cell_list $scan_cell_list 3] * 70/100}]                                    
                                } else {
                                    # penalize it depending on the slack
                                    lset cell_list $scan_cell_list 3 [expr {[lindex $cell_list $scan_cell_list 3] * -1 * [lindex $cell_list $scan_cell_list 6]/$min_slack}]
                                }
                            }
                            set scan_cell_list [expr {$scan_cell_list + 1}]
                        }
                        set explore_path [expr {$explore_path + 1}]
                    }
                    set index [expr {$index + 1}]
                }
                set cell_list [lsort -decreasing -real -index 3 [lrange $cell_list 0 [expr {[llength $cell_list] - 1}]]]
            } else {
                set cell_fixed [lrange $cell_list 0 [expr {($to_change - $step*2)-1}]]
                # Sort the part not fixed
                set sorted_part [lsort -decreasing -real -index 3 [lrange $cell_list [expr {$to_change - $step*2}] [expr {[llength $cell_list]-1}]]]
                set index 0
                while {$index < [llength $sorted_part]} {
                    set sorted_part [lsort -decreasing -real -index 3 [lrange $cell_list [expr {$to_change - $step*2}] [expr {[llength $cell_list]-1}]]]
                    set curr_cell [lindex $sorted_part $index 4]
                    set explore_path 0
                    set path_list [lindex $sorted_part $index 5]
                    while {$explore_path < [llength $path_list]} {
                        set scan_cell_list 0
                        while {$scan_cell_list < [llength $sorted_part]} {
                            if {[lindex $sorted_part $scan_cell_list 4] == [lindex $path_list $explore_path]} {  
                                if {[lindex $sorted_part $scan_cell_list 6] >= 0} {
                                    # penalize it a lot
                                    lset sorted_part $scan_cell_list 3 [expr {[lindex $sorted_part $scan_cell_list 3] * 70/100}]                                  
                                } else {
                                    # penalize it depending on the slack
                                    lset sorted_part $scan_cell_list 3 [expr {[lindex $sorted_part $scan_cell_list 3] * -1 * [lindex $sorted_part $scan_cell_list 6]/$min_slack}]
                                }
                                # decrease the weight (uniformly for all cells in the worst slack path of the currently selected cell)
                            }
                            set scan_cell_list [expr {$scan_cell_list + 1}]
                        }
                        set explore_path [expr {$explore_path + 1}]
                    }
                    set index [expr {$index + 1}]
                }
                set sorted_part [lsort -decreasing -real -index 3 [lrange $cell_list [expr {$to_change - $step*2}] [expr {[llength $cell_list]-1}]]]
                set cell_list [concat $cell_fixed $sorted_part]
            }
        }
        # Sort if you are in the first iteration
        if {$prev_outcome == 0} {
            set index 0
            while {$index < $to_change} {
                set curr_cell [lindex $cell_list $index 4]
                set explore_path 0
                set path_list [lindex $cell_list $index 5]
                while {$explore_path < [llength [lindex $cell_list $index 5]]} {
                    set scan_cell_list $index
                    while {$scan_cell_list < [llength $cell_list]} {
                        if {[lindex $cell_list $scan_cell_list 4] == [lindex $path_list $explore_path]} {  
                            # decrease the weight (uniformly for all cells in the worst slack path of the currently selected cell)
                            if {[lindex $cell_list $scan_cell_list 6] >= 0} {
                                # penalize it a lot
                                lset cell_list $scan_cell_list 3 [expr {[lindex $cell_list $scan_cell_list 3] * 70/100}]                                    
                            } else {
                                # penalize it depending on the slack
                                lset cell_list $scan_cell_list 3 [expr {[lindex $cell_list $scan_cell_list 3] * -1 * [lindex $cell_list $scan_cell_list 6]/$min_slack }]
                            }
                        }
                        set scan_cell_list [expr {$scan_cell_list + 1}]
                    }
                    set explore_path [expr {$explore_path + 1}]
                }
                set index [expr {$index + 1}]
            }
            set cell_list [lsort -decreasing -real -index 3 [lrange $cell_list 0 [expr {[llength $cell_list] - 1}]]]
        }
        # emulation of do while
        while { 0 < 1 } {
            # CHECK FOR MAX ITER
            if {$step == 0} {
                break
            }
            # change the cells
            if {$prev_outcome == 0} { 
                # Initial state
                cells_swapping [lrange [lmap x $cell_list {lindex $x 4}] 0 [expr {$to_change-1}]] LVT
            } elseif {$prev_outcome == 1} {
                # Constraints met in previous iteration
                cells_swapping [lrange [lmap x $cell_list {lindex $x 4}] [expr {$to_change-1}] [expr {($to_change + $step * 2)-1}]] HVT
            } elseif {$prev_outcome == 2} {
                # Constraints NOT met in previous iteration
                cells_swapping [lrange [lmap x $cell_list {lindex $x 4}] 0 [expr {$to_change-1}]] LVT
            }
            set iter [expr {$iter + 1}]
            # check if constraints met
            set ::current_FECs {}
            set ::cells_in_paths {}
            set ::worst_slacks {}
            if { [constraints_met $slackThreshold $maxFanoutEndpointCost $cell_list] == 1} {
                # MET
                set leakage_final [get_attribute [current_design] leakage_power]
                set normalized_power [expr ($leakage_final / $leakage_initial)]
                # set the working_found flag
                set working_found 1

                # reset lists of best configuration
                set lvt_cells {}
                set hvt_cells {}
                # save this configuration
                foreach elem $cell_list {
                    set cell [get_cells [lindex $elem 4]]
                    set original_vt [get_attribute $cell lib_cell.threshold_voltage_group]
                    if {$original_vt == "LVT"} {
                        lappend lvt_cells [lindex $elem 4]
                    } elseif {$original_vt == "HVT"} {
                        lappend hvt_cells [lindex $elem 4]
                    }
                }
                # reduce the number of cells to be changed
                set to_change [expr {$to_change - $step}]
                # update odd_divisions
                if {[expr {$to_change % 2}] == 1 } {
                    set odd_divisions [expr {$odd_divisions + 1}]
                }
                # change the step
                set step [expr {$step / 2}]
                # set the outcome to 1
                set prev_outcome 1
            } else {
                set to_change [expr {$to_change + $step}]
                # update odd_divisions
                if {[expr {$to_change % 2}] == 1 } {
                    set odd_divisions [expr {$odd_divisions + 1}]
                }
                # update the step, you are moving exactly as you would move if the constraints were met
                set step [expr {$step / 2}]
                # set the outcome to 2
                set prev_outcome 2
                # back to step 1
                break    
            }
        }
    }
    set step 1
    set has_entered_adjustment_loop 0
    # to complete the exploration of the array in case no valid solution is found through the search
    while {$to_change <= [llength $cell_list] && $working_found == 0} {
        set has_entered_adjustment_loop 1
        cells_swapping [lrange [lmap x $cell_list {lindex $x 4}] 0 [expr {$to_change-1}]] LVT
        # if constraints are met break and return
        if { [constraints_met $slackThreshold $maxFanoutEndpointCost $cell_list] == 1} {
            return 1
        }
        set to_change [expr {$to_change + $step}]
    }
    if {$has_entered_adjustment_loop == 1} {
        return 0
    }
    
    # REVERT TO THE BEST WORKING SOLUTION (you need it if you want to go on with the binary search on the left LVT array)
    if { $working_found == 1 } {
        # ripristinate best solution local to the current execution
        cells_swapping $lvt_cells LVT
        cells_swapping $hvt_cells HVT
        return 1
    }
    return 0
    
}