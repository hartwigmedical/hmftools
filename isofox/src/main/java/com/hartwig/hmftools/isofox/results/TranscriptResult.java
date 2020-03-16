package com.hartwig.hmftools.isofox.results;

import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.GeneCollection.TRANS_COUNT;
import static com.hartwig.hmftools.isofox.common.GeneCollection.UNIQUE_TRANS_COUNT;
import static com.hartwig.hmftools.isofox.common.RegionReadData.findExonRegion;
import static com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesGenerator.FL_FREQUENCY;
import static com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesGenerator.FL_LENGTH;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStrList;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.isofox.common.FragmentMatchType;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesGenerator;

import org.immutables.value.Value;

@Value.Immutable
public abstract class TranscriptResult
{
    public abstract TranscriptData trans();

    public abstract int exonsFound();
    public abstract int shortSupportingFragments();
    public abstract int shortUniqueFragments();
    public abstract int longSupportingFragments();
    public abstract int longUniqueFragments();
    public abstract int spliceJunctionsSupported();
    public abstract int spliceJunctionFragments();
    public abstract int spliceJunctionUniqueFragments();
    public abstract int uniqueSpliceJunctions();
    public abstract int uniqueSpliceJunctionsSupported();
    public abstract int exonicBases();
    public abstract int exonicBaseCoverage();
    public abstract int uniqueBases();
    public abstract int uniqueBaseCoverage();
    public abstract double uniqueBaseAvgDepth();

    public abstract double effectiveLength();
    public abstract double fitAllocation();

    public static TranscriptResult createTranscriptResults(
            final GeneCollection geneCollection, final GeneReadData geneReadData, final TranscriptData transData,
            final List<int[]> expRateFragmentLengths)
    {
        double expRateAllocation = geneCollection.getFitAllocation(transData.TransName);

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

            final RegionReadData exonReadData = findExonRegion(geneReadData.getExonRegions(), exon.ExonStart, exon.ExonEnd);
            if(exonReadData == null)
                continue;

            int exonCoverage = exonReadData.baseCoverage(1);
            exonicBaseCoverage += exonCoverage;

            if(exonCoverage > 0)
                ++exonsFound;

            if(exonReadData.getTransExonRefs().size() == 1)
            {
                uniqueExonicBases += exonReadData.uniqueBaseCount();
                uniqueExonicBaseCoverage += exonReadData.uniqueBaseCoverage(1);
                uniqueExonicBaseTotalDepth += exonReadData.uniqueBaseTotalDepth();
            }

            exonicBases += exon.ExonEnd - exon.ExonStart + 1;

            if(i > 0)
            {
                int[] sjReads = exonReadData.getTranscriptJunctionMatchCount(transData.TransId, SE_START);

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

        int[][] supportingFragments = geneCollection.getTranscriptReadCount(transData.TransId);

        double effectiveLength = calcEffectiveLength(exonicBases, expRateFragmentLengths);

        if(effectiveLength < 0)
        {
            final List<String> fragDataStr = expRateFragmentLengths.stream()
                    .map(x -> String.format("%d=%d", x[FL_LENGTH], x[FL_FREQUENCY])).collect(Collectors.toList());

            ISF_LOGGER.error("trans({}) exonicBases({}) neg calculated effective length({}) from dist({})",
                    transData.TransName, exonicBases, effectiveLength, appendStrList(fragDataStr, ';'));
        }

        TranscriptResult results = ImmutableTranscriptResult.builder()
                .trans(transData)
                .exonsFound(exonsFound)
                .spliceJunctionsSupported(spliceJunctionsSupported)
                .uniqueSpliceJunctions(uniqueSpliceJunctions)
                .uniqueSpliceJunctionsSupported(uniqueSpliceJunctionsSupported)
                .spliceJunctionFragments(supportingFragments[FragmentMatchType.typeAsInt(FragmentMatchType.SPLICED)][TRANS_COUNT])
                .spliceJunctionUniqueFragments(supportingFragments[FragmentMatchType.typeAsInt(FragmentMatchType.SPLICED)][UNIQUE_TRANS_COUNT])
                .shortSupportingFragments(supportingFragments[FragmentMatchType.typeAsInt(FragmentMatchType.SHORT)][TRANS_COUNT])
                .shortUniqueFragments(supportingFragments[FragmentMatchType.typeAsInt(FragmentMatchType.SHORT)][UNIQUE_TRANS_COUNT])
                .longSupportingFragments(supportingFragments[FragmentMatchType.typeAsInt(FragmentMatchType.LONG)][TRANS_COUNT])
                .longUniqueFragments(supportingFragments[FragmentMatchType.typeAsInt(FragmentMatchType.LONG)][UNIQUE_TRANS_COUNT])
                .exonicBases(exonicBases)
                .exonicBaseCoverage(exonicBaseCoverage)
                .uniqueBases(uniqueExonicBases)
                .uniqueBaseCoverage(uniqueExonicBaseCoverage)
                .uniqueBaseAvgDepth(uniqueBaseAvgDepth)
                .effectiveLength(effectiveLength)
                .fitAllocation(expRateAllocation)
                .build();

        return results;
    }

    public static double calcEffectiveLength(int transLength, final List<int[]> fragmentLengthData)
    {
        if(fragmentLengthData.isEmpty())
            return transLength;

        long flFrequencyTotal = 0;
        long flBasesTotal = 0;

        for(final int[] flData : fragmentLengthData)
        {
            int fragLength = flData[FL_LENGTH];
            int fragFrequency = flData[FL_FREQUENCY];

            if(fragLength >= transLength)
                continue;

            long possibleBases = transLength - fragLength;
            flFrequencyTotal += fragFrequency;
            flBasesTotal += possibleBases * fragFrequency;
        }

        if(flFrequencyTotal == 0)
            return 0;

        return flBasesTotal / (double)flFrequencyTotal;
    }

    private static boolean isSpliceJunctionUnique(final String transId, final List<TranscriptData> transDataList, long exonEnd, long exonStart)
    {
        for(TranscriptData transData : transDataList)
        {
            if (transData.TransName.equals(transId))
                continue;

            for (int i = 1; i < transData.exons().size(); ++i)
            {
                if(transData.exons().get(i-1).ExonEnd == exonEnd && transData.exons().get(i).ExonStart == exonStart)
                    return false;
            }
        }

        return true;
    }
}
