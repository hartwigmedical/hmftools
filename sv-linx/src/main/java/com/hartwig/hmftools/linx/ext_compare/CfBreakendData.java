package com.hartwig.hmftools.linx.ext_compare;

import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.apache.commons.compress.utils.Lists;

public class CfBreakendData
{
    public final int ChainId;
    public final int BreakpointId;
    public final int RearrangeId;
    public String Chromosome; // not final since can be set to Nan for X or Y and then corrected
    public final int Position;
    public final byte Orientation;
    public final boolean[] ChromosomeEndLoss;
    public final int DbBreakendId;
    public final List<Integer> AdjacentBreakendIds;

    public static final int NO_ID = -1;

    public CfBreakendData(int chainId, int breakpointId, int rearrangeId, final String chromosome, int position,
            byte orientation, final boolean[] chromosomeEndLoss,
            int dbBreakendId, final List<Integer> adjacentBreakendIds)
    {
        ChainId = chainId;
        BreakpointId = breakpointId;
        RearrangeId = rearrangeId;
        Chromosome = chromosome;
        Position = position;
        Orientation = orientation;
        ChromosomeEndLoss = chromosomeEndLoss;
        DbBreakendId = dbBreakendId;
        AdjacentBreakendIds = adjacentBreakendIds;
    }

    public static final String CF_DATA_DELIMITER = "\t\\s*";

    public static CfBreakendData fromData(final String data, final Map<String,Integer> fieldIndexMap)
    {
        // Sample	 Chain	 Breakpoint	 Rearrangement	 Chromosome:position	 Strand
        // Potential chromosomal end loss (left)	 Potential chromosomal end loss (right)	 Deletion bridge partner breakpoint
        // Adjacent breakpoint(s)	 Site annotation
        //CPCT02010452T	 3	 196	 360203	 14:49979403	 1	 0	 0	 195	 |195|197|	 na
        //CPCT02010452T	 3	 197	 360204	 14:49983004	 0	 0	 0	 198	 |196|198|	 na

        String[] items = data.split(CF_DATA_DELIMITER);

        final int chainId = Integer.parseInt(items[fieldIndexMap.get("Chain")]);
        final int breakpointId = Integer.parseInt(items[fieldIndexMap.get("Breakpoint")]);
        final int rearrangeId = Integer.parseInt(items[fieldIndexMap.get("Rearrangement")]);

        final String breakendData = items[fieldIndexMap.get("Chromosome:position")];
        final String chromosome = breakendData.split(":")[0];
        final int position = Integer.parseInt(breakendData.split(":")[1]);

        final byte orientation = items[fieldIndexMap.get("Strand")].equals("0") ? POS_STRAND : NEG_STRAND;

        final boolean[] chromosomeEndLoss = new boolean[] {
                items[fieldIndexMap.get("Potential chromosomal end loss (left)")].equals("1"),
                items[fieldIndexMap.get("Potential chromosomal end loss (right)")].equals("1") };

        final int dbBreakendId = items.length == fieldIndexMap.size() ?
                Integer.parseInt(items[fieldIndexMap.get("Deletion bridge partner breakpoint")]) : NO_ID;

        int adjBndIndex = fieldIndexMap.get("Adjacent breakpoint(s)");
        final String adjacentBreakendData = items.length == fieldIndexMap.size() ? items[adjBndIndex] : items[adjBndIndex - 1];
        int abLength = adjacentBreakendData.length();

        final List<Integer> adjacentBreakendIds = abLength > 2 ?
                Arrays.stream(adjacentBreakendData.substring(0, abLength - 1).substring(1).split("\\|"))
                        .map(x -> Integer.parseInt(x)).collect(Collectors.toList()) : Lists.newArrayList();

        return new CfBreakendData(chainId, breakpointId, rearrangeId, chromosome, position, orientation, chromosomeEndLoss, dbBreakendId, adjacentBreakendIds);
    }


    // Sample	 Chain	 Breakpoint	 Rearrangement	 Chromosome:position	 Strand
    // Potential chromosomal end loss (left)	 Potential chromosomal end loss (right)	 Deletion bridge partner breakpoint	 Adjacent breakpoint(s)	 Site annotation


}
