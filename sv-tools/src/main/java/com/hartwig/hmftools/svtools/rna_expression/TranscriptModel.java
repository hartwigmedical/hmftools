package com.hartwig.hmftools.svtools.rna_expression;

import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TC_LONG;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TC_SHORT;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TC_SPLICE;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TRANS_COUNT;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.UNIQUE_TRANS_COUNT;
import static com.hartwig.hmftools.svtools.rna_expression.RegionReadData.extractTransId;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class TranscriptModel
{
    private static final Logger LOGGER = LogManager.getLogger(TranscriptModel.class);

    public static TranscriptResults calculateTranscriptResults(final GeneReadData geneReadData, final TranscriptData transData)
    {
        int exonsFound = 0;

        int spliceJunctionsSupported = 0;
        int exonicBases = 0;
        int exonicBaseCoverage = 0;

        int uniqueExonicBases = 0;
        int uniqueExonicBaseCoverage = 0;
        int uniqueExonicBaseTotalDepth = 0;

        int uniqueSpliceJunctions = 0;
        int uniqueSpliceJunctionsSupported = 0;

        /* Criteria for transcript selection
        - all exon junctions covered
        - unique exon junctions
        - split reads skipping exons
        - unique exon reads (but could cover introns as well
         */

        final List<ExonData> exons = transData.exons();

        for(int i = 0; i < exons.size(); ++i)
        {
            ExonData exon = exons.get(i);

            final RegionReadData exonReadData = geneReadData.findExonRegion(exon.ExonStart, exon.ExonEnd);
            if(exonReadData == null)
                continue;

            int exonCoverage = exonReadData.baseCoverage(1);
            exonicBaseCoverage += exonCoverage;

            if(exonCoverage > 0)
                ++exonsFound;

            if(exonReadData.getRefRegions().size() == 1)
            {
                uniqueExonicBases += exonReadData.uniqueBaseCount();
                uniqueExonicBaseCoverage += exonReadData.uniqueBaseCoverage(1);
                uniqueExonicBaseTotalDepth += exonReadData.uniqueBaseTotalDepth();
            }

            exonicBases += exon.ExonEnd - exon.ExonStart + 1;

            if(i > 0)
            {
                int[] sjReads = exonReadData.getTranscriptJunctionMatchCount(transData.TransName, SE_START);

                final ExonData prevExon = exons.get(i - 1);
                boolean sjUnique = isSpliceJunctionUnique(transData.TransName, geneReadData.getTranscripts(), prevExon.ExonEnd, exon.ExonStart);

                if(sjUnique)
                    ++uniqueSpliceJunctions;

                if(sjReads[TRANS_COUNT] > 0)
                {
                    ++spliceJunctionsSupported;

                    if(sjUnique)
                        ++uniqueSpliceJunctionsSupported;
                }
            }
        }

        double uniqueBaseAvgDepth = uniqueExonicBases > 0 ? uniqueExonicBaseTotalDepth / (double)uniqueExonicBases : 0;

        int[][] supportingFragments = geneReadData.getTranscriptReadCount(transData.TransName);

        TranscriptResults results = ImmutableTranscriptResults.builder()
                .trans(transData)
                .exonsFound(exonsFound)
                .spliceJunctionsSupported(spliceJunctionsSupported)
                .uniqueSpliceJunctions(uniqueSpliceJunctions)
                .uniqueSpliceJunctionsSupported(uniqueSpliceJunctionsSupported)
                .spliceJunctionFragments(supportingFragments[TC_SPLICE][TRANS_COUNT])
                .spliceJunctionUniqueFragments(supportingFragments[TC_SPLICE][UNIQUE_TRANS_COUNT])
                .shortSupportingFragments(supportingFragments[TC_SHORT][TRANS_COUNT])
                .shortUniqueFragments(supportingFragments[TC_SHORT][UNIQUE_TRANS_COUNT])
                .longSupportingFragments(supportingFragments[TC_LONG][TRANS_COUNT])
                .longUniqueFragments(supportingFragments[TC_LONG][UNIQUE_TRANS_COUNT])
                .exonicBases(exonicBases)
                .exonicBaseCoverage(exonicBaseCoverage)
                .uniqueBases(uniqueExonicBases)
                .uniqueBaseCoverage(uniqueExonicBaseCoverage)
                .uniqueBaseAvgDepth(uniqueBaseAvgDepth)
                .build();

        return results;
    }

    private static boolean isSpliceJunctionUnique(final String transId, final List<TranscriptData> transDataList, long exonEnd, long exonStart)
    {
        for(TranscriptData transData : transDataList)
        {
            if (transData.TransName.equals(transId))
                continue;

            for (int i = 1; i < transData.exons().size(); ++i)
            {
                if(transData.exons().get(0).ExonEnd == exonEnd && transData.exons().get(1).ExonStart == exonStart)
                    return false;
            }
        }

        return true;
    }


}
