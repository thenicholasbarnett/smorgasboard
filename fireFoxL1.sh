#!/bin/bash

run=387878
#single/double egs span 179 to 185
bitVals=(185 179 181 182 183 136 140 143 145 147 272 276 270)

bitURLStr=""

pos=0
for i in ${bitVals[@]}
do
    #URL STR
    urlStr="https://cmsoms.cern.ch/cms/triggers/l1_algo?cms_run=$run&cms_l1_bit=$i"
    #URL STR w/ lumi restriction
    
#    if [[ $pos -eq 0 ]]
#    then
#	firefox --new-window $urlStr
#    else
#	firefox --new-tab $urlStr
#    fi

    bitURLStr=$bitURLStr"$i":16,
    
    pos=$((pos + 1))
done

bitURLStr=${bitURLStr%,}

firefox --new-tab "https://cmsoms.cern.ch/cms/triggers/l1_rates?cms_run=$run&props.20639_19388.selectedCells=Physics:2&props.19391_19388.selectedCells=L1A%20physics:2&props.19392_19388.selectedCells=$bitURLStr"

#LUMIVERSION
#firefox --new-tab "https://cmsoms.cern.ch/cms/triggers/l1_rates?cms_run=$run&props.20639_19388.selectedCells=Physics:2&props.19391_19388.selectedCells=L1A%20physics:2&props.19392_19388.selectedCells=$bitURLStr"&cms_ls_from=80&cms_ls_to=90
