package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.LILAC_FILE_HLA_Y_COVERAGE;
import static com.hartwig.hmftools.lilac.LilacConstants.LILAC_FILE_HLA_Y_FRAGMENTS;
import static com.hartwig.hmftools.lilac.ReferenceData.getAminoAcidExonBoundaries;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.HLA_Y;
import static com.hartwig.hmftools.lilac.fragment.FragmentSource.REFERENCE;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_A;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.fragment.FragmentSource;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaGene;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.seq.SequenceMatchType;

public class HlaYCoverage
{
    private final List<HlaSequenceLoci> mHlaYSequences;
    private final Set<Integer> mAminoAcidHetLoci;
    private final List<FragmentAlleles> mRemovedRefFragAlleles;

    private boolean mExceedsThreshold;

    private final Map<FragmentSource, Map<HlaAllele, int[]>> mSourceAlleleFragmentCounts;
    private final Map<FragmentSource, int[]> mSourceMiscCounts;

    private BufferedWriter mWriter;
    private final LilacConfig mConfig;

    private static final int UNIQUE = 0;
    private static final int SHARED_HLAY = 1; // shared with 2+ of the HLA-Y alleles
    private static final int SHARED = 2; // shared with one or more of the solution alleles
    private static final int MAX = SHARED + 1;

    private static final int Y0101_X = 0;
    private static final int EXON_3 = 1;

    public HlaYCoverage(
            final List<HlaSequenceLoci> hlaYSequences, final Map<HlaGene, Map<Integer, Set<String>>> geneAminoAcidHetLociMap,
            final LilacConfig config)
    {
        mHlaYSequences = hlaYSequences;
        mAminoAcidHetLoci = Sets.newHashSet(geneAminoAcidHetLociMap.get(HLA_A).keySet());
        mConfig = config;

        mSourceAlleleFragmentCounts = Maps.newHashMap();
        mSourceMiscCounts = Maps.newHashMap();
        mRemovedRefFragAlleles = Lists.newArrayList();
        mExceedsThreshold = false;
    }

    public boolean exceedsThreshold() { return mExceedsThreshold; }

    public HlaAllele getSelectedAllele()
    {
        if(!mExceedsThreshold)
            return null;

        Map<HlaAllele, int[]> alleleCounts = mSourceAlleleFragmentCounts.get(REFERENCE);

        if(alleleCounts == null)
            return null;

        HlaAllele maxAllele = null;
        int maxUniqueCount = 0;

        for(Map.Entry<HlaAllele, int[]> entry : alleleCounts.entrySet())
        {
            int uniqueCount = entry.getValue()[UNIQUE];

            if(uniqueCount > 0 && uniqueCount > maxUniqueCount)
            {
                maxUniqueCount = uniqueCount;
                maxAllele = entry.getKey();
            }
        }

        return maxAllele;
    }

    public void updateAminoAcidLoci(final Map<HlaGene, Map<Integer, Set<String>>> geneAminoAcidHetLociMap)
    {
        mAminoAcidHetLoci.clear();
        mAminoAcidHetLoci.addAll(geneAminoAcidHetLociMap.get(HLA_A).keySet());
    }

    public void checkThreshold(final List<FragmentAlleles> fragAlleles, final Collection<Fragment> fragments)
    {
        // test for presence of HLA-Y and strip out from consideration any fragment mapping to it
        int uniqueHlaY = 0;

        List<FragmentAlleles> matchedFragmentAlleles = Lists.newArrayList();

        // ignore fragments which don't contain any heterozygous locations
        // only test heterozygous locations in A since HLA-Y matches its exon boundaries
        for(Fragment fragment : fragments)
        {
            List<Integer> fragAminoAcidLoci = fragment.aminoAcidsByLoci().keySet().stream().filter(mAminoAcidHetLoci::contains).toList();

            if(fragAminoAcidLoci.isEmpty())
                continue;

            Set<Integer> fragNucleotideLoci = fragment.nucleotidesByLoci().keySet();

            boolean matchesY = false;
            FragmentAlleles matchedFrag = null;

            for(HlaSequenceLoci sequence : mHlaYSequences)
            {
                String fragNucleotides = fragment.nucleotides(fragNucleotideLoci);

                SequenceMatchType matchType = sequence.determineMatchType(fragNucleotides, fragNucleotideLoci);
                if(matchType == SequenceMatchType.FULL)
                {
                    matchesY = true;

                    matchedFrag = fragAlleles.stream()
                            .filter(x -> x.getFragment().id().equals(fragment.id())).findFirst().orElse(null);

                    break;
                }
            }

            if(!matchesY)
                continue;

            if(matchedFrag == null)
            {
                ++uniqueHlaY;
                fragment.setScope(HLA_Y, true);
            }
            else
            {
                matchedFragmentAlleles.add(matchedFrag);
            }
        }

        int totalHlaYFrags = uniqueHlaY + matchedFragmentAlleles.size();
        double threshold = fragments.size() * mConfig.HlaYPercentThreshold;
        mExceedsThreshold = uniqueHlaY >= threshold;

        if(totalHlaYFrags > 0)
        {
            LL_LOGGER.info("HLA-Y fragments({} unique={}) shared={}) {}",
                    totalHlaYFrags, uniqueHlaY, matchedFragmentAlleles.size(),
                    mExceedsThreshold ? "above threshold" : "below threshold");

            if(mExceedsThreshold)
            {
                matchedFragmentAlleles.forEach(fragAlleles::remove);

                mRemovedRefFragAlleles.addAll(matchedFragmentAlleles); // remember for later counting purposes
                mRemovedRefFragAlleles.forEach(x -> x.getFragment().setScope(HLA_Y, true));
            }
        }
    }

    public void assignReferenceFragments(
            final Collection<HlaAllele> alleles, final Iterable<FragmentAlleles> fragAlleles, final Iterable<Fragment> fragments)
    {
        // put back those fragments previously excluded due to being matched to HLA-Y
        List<FragmentAlleles> allFragAlleles = Lists.newArrayList(fragAlleles);
        allFragAlleles.addAll(mRemovedRefFragAlleles);

        assignFragments(alleles, allFragAlleles, fragments, REFERENCE);
    }

    public void assignFragments(
            final Collection<HlaAllele> alleles, final List<FragmentAlleles> fragAlleles, final Iterable<Fragment> fragments,
            final FragmentSource source)
    {
        Map<HlaAllele, int[]> sourceAlleleCounts = Maps.newHashMap();
        mSourceAlleleFragmentCounts.put(source, sourceAlleleCounts);

        int[] miscCounts = new int[2];
        mSourceMiscCounts.put(source, miscCounts);

        List<FragmentAlleles> matchedFragmentAlleles = Lists.newArrayList();

        int uniqueHlaY = 0;

        for(Fragment fragment : fragments)
        {
            List<Integer> fragAminoAcidLoci = fragment.aminoAcidsByLoci().keySet().stream()
                    .filter(mAminoAcidHetLoci::contains)
                    .collect(Collectors.toList());

            if(fragAminoAcidLoci.isEmpty())
                continue;

            Set<Integer> fragNucleotideLoci = fragment.nucleotidesByLoci().keySet();

            List<HlaAllele> matchedAlleles = Lists.newArrayList();

            for(HlaSequenceLoci sequence : mHlaYSequences)
            {
                String fragNucleotides = fragment.nucleotides(fragNucleotideLoci);

                SequenceMatchType matchType = sequence.determineMatchType(fragNucleotides, fragNucleotideLoci);

                if(matchType == SequenceMatchType.FULL)
                    matchedAlleles.add(sequence.Allele);
            }

            if(matchedAlleles.isEmpty())
                continue;

            // this fragment supports 1 or more HLA-Y alleles, now check if it is also used in one of the solution alleles
            // and capture of key metrics
            updateMiscCounts(fragment, miscCounts);

            FragmentAlleles matchedFragment = fragAlleles.stream()
                    .filter(x -> x.getFragment().id().equals(fragment.id()))
                    .filter(x -> fragmentSupportsSolution(x, alleles))
                    .findFirst().orElse(null);

            if(matchedFragment != null)
                matchedFragmentAlleles.add(matchedFragment);
            else
                ++uniqueHlaY;

            if(LL_LOGGER.isDebugEnabled())
            {
                writeFragmentInfo(fragment, fragAminoAcidLoci, matchedAlleles, matchedFragment != null);
            }

            for(HlaAllele allele : matchedAlleles)
            {
                int[] alleleCounts = sourceAlleleCounts.get(allele);

                if(alleleCounts == null)
                {
                    alleleCounts = new int[MAX];
                    sourceAlleleCounts.put(allele, alleleCounts);
                }

                if(matchedFragment != null)
                    ++alleleCounts[SHARED];
                else if(matchedAlleles.size() > 1)
                    ++alleleCounts[SHARED_HLAY];
                else
                    ++alleleCounts[UNIQUE];
            }
        }

        if(source != REFERENCE)
        {
            LL_LOGGER.info("HLA-Y src({}) fragments({} unique={}) shared={})",
                    source, uniqueHlaY + matchedFragmentAlleles.size(), uniqueHlaY, matchedFragmentAlleles.size());

            if(mExceedsThreshold)
                matchedFragmentAlleles.forEach(fragAlleles::remove);
        }
    }

    private static boolean fragmentSupportsSolution(final FragmentAlleles fragmentAlleles, final Collection<HlaAllele> alleles)
    {
        return alleles.stream().anyMatch(fragmentAlleles::contains);
    }

    private static final int Y0101_X_LOCUS = 227;

    private static void updateMiscCounts(final Fragment fragment, final int[] miscCounts)
    {
        if(fragment.aminoAcidsByLoci().containsKey(Y0101_X_LOCUS) && fragment.aminoAcid(Y0101_X_LOCUS).equals("X"))
            ++miscCounts[Y0101_X];

        List<Integer> hlaAExonBoundaries = getAminoAcidExonBoundaries(HLA_A);
        int exon3Start = hlaAExonBoundaries.get(1) + 1;
        int exon3End = hlaAExonBoundaries.get(2);

        if(fragment.minAminoAcidLocus() >= exon3Start && fragment.maxAminoAcidLocus() <= exon3End)
            ++miscCounts[EXON_3];
    }

    public void writeAlleleCounts(final String sampleId)
    {
        if(mConfig.OutputDir.isEmpty() || !mExceedsThreshold)
            return;

        for(Map.Entry<FragmentSource, int[]> entry : mSourceMiscCounts.entrySet())
        {
            FragmentSource source = entry.getKey();
            final int[] miscCounts = entry.getValue();

            String outcome = (miscCounts[Y0101_X] == 0 && miscCounts[EXON_3] > 0) ? "NOVEL" : "REF";

            LL_LOGGER.info("HLA-Y_MISC_COUNTS:{},{},{},{},{}",
                    sampleId, source, miscCounts[Y0101_X], miscCounts[EXON_3], outcome);
        }

        String fileName = mConfig.formFileId(LILAC_FILE_HLA_Y_COVERAGE, false);
        try
        {
            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write("Source\tAllele\tTotal\tShared\tHlaYShared\tUnique");
            writer.newLine();

            for(Map.Entry<FragmentSource, Map<HlaAllele, int[]>> sourceEntry : mSourceAlleleFragmentCounts.entrySet())
            {
                FragmentSource source = sourceEntry.getKey();
                for(Map.Entry<HlaAllele, int[]> alleleEntry : sourceEntry.getValue().entrySet())
                {
                    final HlaAllele allele = alleleEntry.getKey();
                    final int[] counts = alleleEntry.getValue();

                    int total = counts[SHARED] + counts[SHARED_HLAY] + counts[UNIQUE];

                    writer.write(String.format("%s\t%s\t%d\t%d\t%d\t%d",
                            source, allele, total, counts[SHARED], counts[SHARED_HLAY], counts[UNIQUE]));

                    writer.newLine();
                }
            }

            writer.close();

            closeBufferedWriter(mWriter);
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to write HLA-Y coverage file({}): {}", fileName, e.toString());
        }
    }

    private void writeFragmentInfo(final Fragment fragment, final Iterable<Integer> fragAminoAcidLoci,
            final Iterable<HlaAllele> hlayAlleles, boolean isShared)
    {
        if(mConfig.OutputDir.isEmpty() || !mExceedsThreshold)
            return;

        try
        {
            if(mWriter == null)
            {
                String fileName = mConfig.formFileId(LILAC_FILE_HLA_Y_FRAGMENTS, false);
                mWriter = createBufferedWriter(fileName, false);

                mWriter.write("ReadId\tReadInfo\tAlleles\tShared\tLoci");
                mWriter.newLine();
            }

            StringJoiner sjAlleles = new StringJoiner(ITEM_DELIM);
            hlayAlleles.forEach(x -> sjAlleles.add(x.toString()));

            StringJoiner sjLoci = new StringJoiner(ITEM_DELIM);
            fragAminoAcidLoci.forEach(x -> sjLoci.add(x.toString()));

            mWriter.write(String.format("%s\t%s\t%s\t%s\t%s",
                    fragment.id(), fragment.readInfo(), sjAlleles, isShared, sjLoci));
            mWriter.newLine();

        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to write HLA-Y fragments file: {}", e.toString());
        }
    }
}
