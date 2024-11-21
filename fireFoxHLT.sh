#!/bin/bash

#PRIMARY URLSTR
urlStr="https://cmsoms.cern.ch/cms/triggers/hlt_trigger_rates?props.21024_21019.selectedCells=Physics:2&props.21026_21019.selectedCells=Physics*:2&props.21021_21019.selectedCells=INSERTTRIG&cms_run=INSERTRUN"

#LUMI-restricted URL STR
#urlStr="https://cmsoms.cern.ch/cms/triggers/hlt_trigger_rates?props.21024_21019.selectedCells=Physics:2&props.21026_21019.selectedCells=Physics*:2&props.21021_21019.selectedCells=INSERTTRIG&cms_run=INSERTRUN"

hltStr=""
runStr=387892
triggersJet=(HLT_HIPuAK4CaloJet40Eta5p1_MinBiasHF1AND_ZDC1nOR_v4 HLT_HIPuAK4CaloJet40Eta5p1_MinBiasHF1AND_ZDC2nOR_v6 HLT_HIPuAK4CaloJet40Eta5p1_MinBiasHF1AND_v6 HLT_HIPuAK4CaloJet60Eta5p1_MinBiasHF1AND_ZDC1nOR_v4 HLT_HIPuAK4CaloJet60Eta5p1_MinBiasHF1AND_ZDC2nOR_v6 HLT_HIPuAK4CaloJet60Eta5p1_MinBiasHF1AND_v6 HLT_HIPuAK4CaloJet80Eta5p1_v14 HLT_HIPuAK4CaloJet100Eta5p1_v14 HLT_HIPuAK4CaloJet120Eta2p1_v7 HLT_HIPuAK4CaloJet120Eta5p1_v14 HLT_HIPuAK4CaloJet40Fwd_v7 HLT_HIPuAK4CaloJet60Fwd_v7 HLT_HIPuAK4CaloJet80Fwd_v7 HLT_HIPuAK4CaloJet100Fwd_v7 HLT_HIPuAK4CaloJet120Fwd_v7)

for i in ${triggersJet[@]}
do
    hltStr=$hltStr"$i":2,
done

hltStr=${hltStr%,}

jetURLStr=$(echo $urlStr | sed -e "s@INSERTRUN@$runStr@g")
jetURLStr=$(echo $jetURLStr | sed -e "s@INSERTTRIG@$hltStr@g")

echo $jetURLStr

echo ""

triggersEG=(HLT_HIDoubleEle10GsfMass50_v14 HLT_HIDoubleEle10Gsf_v14 HLT_HIDoubleEle15GsfMass50_v14 HLT_HIDoubleEle15Gsf_v14 HLT_HIDoubleGEDPhoton20_v7 HLT_HIEle10Gsf_v14 HLT_HIEle15Ele10GsfMass50_v14 HLT_HIEle15Ele10Gsf_v14 HLT_HIEle15Gsf_v14 HLT_HIEle20Gsf_v14 HLT_HIEle30Gsf_v14 HLT_HIEle40Gsf_v14 HLT_HIEle50Gsf_v14 HLT_HIGEDPhoton10_EB_v14 HLT_HIGEDPhoton10_v14 HLT_HIGEDPhoton20_EB_v14 HLT_HIGEDPhoton20_v14 HLT_HIGEDPhoton30_EB_v14 HLT_HIGEDPhoton30_v14 HLT_HIGEDPhoton40_EB_v14 HLT_HIGEDPhoton40_v14 HLT_HIGEDPhoton50_EB_v14 HLT_HIGEDPhoton50_v14 HLT_HIGEDPhoton60_EB_v14 HLT_HIGEDPhoton60_v14) 

hltStr=""
for i in ${triggersEG[@]}
do
    hltStr=$hltStr"$i":2,
done

hltStr=${hltStr%,}

egURLStr=$(echo $urlStr | sed -e "s@INSERTRUN@$runStr@g")
egURLStr=$(echo $egURLStr | sed -e "s@INSERTTRIG@$hltStr@g")

echo $egURLStr

open /Applications/Firefox.app -u $egURLStr
open /Applications/Firefox.app -u $jetURLStr
#firefox -u $egURLStr
#firefox -u $jetURLStr
