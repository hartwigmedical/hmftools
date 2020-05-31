package com.hartwig.hmftools.linx.ext_compare;

import java.util.List;

public class CfSummaryData
{
    /*
    - Total breakpoints: 232
- Number of breakpoints in chains: 108
- Percentage of breakpoints assigned to chains: 4.655172e-01
- Number of chains: 5
- Number of interchromosomal chain fusions: 26
- Number of intrachromosomal chain fusions: 82
- Percent of chain fusions that are intrachromosomal: 2.407407e-01
- Max number of interchromosomal fusions in a chain: 24
- Maximum number of chromosomes involved in one chain: 4
- Number of breakpoints in longest chain: 70
- Number of deletion bridges: 16
     */

    public final int TotalBreakpoints;
    public final int BreakendsInChains;
    public final int ChainCount;
    public final int MultiChromosomeFusions;
    public final int SingleChromosomeFusions;
    public final int MaxChromosomeChain;
    public final int MaxChainCount;
    public final int DeletionBridgeCount;

    public CfSummaryData()
    {
        TotalBreakpoints = 0;
        BreakendsInChains = 0;
        ChainCount = 0;
        MultiChromosomeFusions = 0;
        SingleChromosomeFusions = 0;
        MaxChromosomeChain = 0;
        MaxChainCount = 0;
        DeletionBridgeCount = 0;
    }

    public void fromDataFile(final List<String> data)
    {

    }

}
