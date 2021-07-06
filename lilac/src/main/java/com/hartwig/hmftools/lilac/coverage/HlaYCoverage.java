package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.A_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_A;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_Y_FRAGMENT_THRESHOLD;
import static com.hartwig.hmftools.lilac.LilacConstants.ITEM_DELIM;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.HLA_Y;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.seq.SequenceMatchType;

public class HlaYCoverage
{
    private final List<HlaSequenceLoci> mHlaYSequences;
    private final Set<Integer> mAminoAcidHetLoci;
    private final List<FragmentAlleles> mRemovedRefFragAlleles;

    private boolean mExceedsThreshold;
    private final Map<String,Map<HlaAllele,int[]>> mSourceAlleleFragmentCounts;
    private final Map<String,int[]> mSourceMiscCounts;

    private BufferedWriter mWriter;
    private final String mOutputPrefix;

    private static final int UNIQUE = 0;
    private static final int SHARED_HLAY = 1; // shared with 2+ of the HLA-Y alleles
    private static final int SHARED = 2; // shared with one or more of the solution alleles
    private static final int MAX = SHARED + 1;

    private static final int Y0101_X = 0;
    private static final int EXON_3 = 1;

    private static final String REF_SOURCE = "REF";

    public HlaYCoverage(
            final List<HlaSequenceLoci> hlaYSequences, final Map<String,Map<Integer,Set<String>>> geneAminoAcidHetLociMap,
            final String outputPrefix)
    {
        mHlaYSequences = hlaYSequences;
        mAminoAcidHetLoci = Sets.newHashSet(geneAminoAcidHetLociMap.get(GENE_A).keySet());
        mOutputPrefix = outputPrefix;

        mSourceAlleleFragmentCounts = Maps.newHashMap();
        mSourceMiscCounts = Maps.newHashMap();
        mRemovedRefFragAlleles = Lists.newArrayList();
        mExceedsThreshold = false;
    }

    public boolean exceedsThreshold() { return mExceedsThreshold; }

    public void updateAminoAcidLoci(final Map<String,Map<Integer,Set<String>>> geneAminoAcidHetLociMap)
    {
        mAminoAcidHetLoci.clear();
        mAminoAcidHetLoci.addAll(geneAminoAcidHetLociMap.get(GENE_A).keySet());
    }

    public void checkThreshold(final List<FragmentAlleles> fragAlleles, final List<Fragment> fragments)
    {
        // test for presence of HLA-Y and strip out from consideration any fragment mapping to it
        int uniqueHlaY = 0;

        List<FragmentAlleles> matchedFragmentAlleles = Lists.newArrayList();

        // ignore fragments which don't contain any heterozygous locations
        // only test heterozygous locations in A since HLA-Y matches its exon boundaries
        for(Fragment fragment : fragments)
        {
            List<Integer> fragAminoAcidLoci = fragment.getAminoAcidLoci().stream()
                    .filter(x -> mAminoAcidHetLoci.contains(x)).collect(Collectors.toList());

            if(fragAminoAcidLoci.isEmpty())
                continue;

            List<Integer> fragNucleotideLoci = fragment.getNucleotideLoci();

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

        int totalHlaYFrags = uniqueHlaY  + matchedFragmentAlleles.size();
        double threshold = fragments.size() * HLA_Y_FRAGMENT_THRESHOLD;
        mExceedsThreshold = uniqueHlaY >= threshold;

        if(totalHlaYFrags > 0)
        {
            LL_LOGGER.info("HLA-Y fragments({} unique={}) shared={}) {}",
                    totalHlaYFrags, uniqueHlaY, matchedFragmentAlleles.size(),
                    mExceedsThreshold ? "above threshold" : "below threshold");

            if(mExceedsThreshold)
            {
                matchedFragmentAlleles.forEach(x -> fragAlleles.remove(x));

                mRemovedRefFragAlleles.addAll(matchedFragmentAlleles); // remember for later counting purposes
                mRemovedRefFragAlleles.forEach(x -> x.getFragment().setScope(HLA_Y, true));
            }
        }
    }

    public void assignReferenceFragments(
            final List<HlaAllele> alleles, final List<FragmentAlleles> fragAlleles, final List<Fragment> fragments)
    {
        // put back those fragments previously excluded due to being matched to HLA-Y
        List<FragmentAlleles> allFragAlleles = Lists.newArrayList(fragAlleles);
        allFragAlleles.addAll(mRemovedRefFragAlleles);

        assignFragments(alleles, allFragAlleles, fragments, REF_SOURCE);
    }

    public void assignFragments(
            final List<HlaAllele> alleles, final List<FragmentAlleles> fragAlleles, final List<Fragment> fragments, final String source)
    {
        Map<HlaAllele,int[]> sourceAlleleCounts = Maps.newHashMap();
        mSourceAlleleFragmentCounts.put(source, sourceAlleleCounts);

        int[] miscCounts = new int[2];
        mSourceMiscCounts.put(source, miscCounts);

        List<FragmentAlleles> matchedFragmentAlleles = Lists.newArrayList();

        int uniqueHlaY = 0;

        for(Fragment fragment : fragments)
        {
            List<Integer> fragAminoAcidLoci = fragment.getAminoAcidLoci().stream()
                    .filter(x -> mAminoAcidHetLoci.contains(x)).collect(Collectors.toList());

            if(fragAminoAcidLoci.isEmpty())
                continue;

            List<Integer> fragNucleotideLoci = fragment.getNucleotideLoci();

            updateMiscCounts(miscCounts, fragAminoAcidLoci);

            FragmentAlleles matchedFragment = fragAlleles.stream()
                    .filter(x -> x.getFragment().id().equals(fragment.id()))
                    .filter(x -> fragmentSupportsSolution(x, alleles))
                    .findFirst().orElse(null);

            if(matchedFragment != null)
                matchedFragmentAlleles.add(matchedFragment);
            else
                ++uniqueHlaY;

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

            writeFragmentInfo(fragment, fragAminoAcidLoci, matchedAlleles, matchedFragment != null);

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

        if(!source.equals(REF_SOURCE))
        {
            LL_LOGGER.info("HLA-Y src({}) fragments({} unique={}) shared={})",
                    source, uniqueHlaY + matchedFragmentAlleles.size(), uniqueHlaY, matchedFragmentAlleles.size());

            matchedFragmentAlleles.forEach(x -> fragAlleles.remove(x));
        }
    }

    private static boolean fragmentSupportsSolution(final FragmentAlleles fragmentAlleles, final List<HlaAllele> alleles)
    {
        return alleles.stream().anyMatch(x -> fragmentAlleles.contains(x));
    }

    public void writeAlleleCounts()
    {
        if(mOutputPrefix == null || !mExceedsThreshold)
            return;

        for(Map.Entry<String,int[]> entry : mSourceMiscCounts.entrySet())
        {
            final String source = entry.getKey();
            final int[] miscCounts = entry.getValue();

            String outcome = (miscCounts[Y0101_X] == 0 && miscCounts[EXON_3] == 0) ? "NOVEL" : "REF";

            LL_LOGGER.info("HLA-Y_MISC_COUNTS:{},{},{},{}",
                    source, miscCounts[Y0101_X], miscCounts[EXON_3], outcome);
        }

        String fileName = mOutputPrefix + ".hlay.coverage.csv";

        try
        {
            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write("Source,Allele,Total,Shared,HlaYShared,Unique");
            writer.newLine();

            for(Map.Entry<String,Map<HlaAllele,int[]>> sourceEntry : mSourceAlleleFragmentCounts.entrySet())
            {
                final String source = sourceEntry.getKey();
                for(Map.Entry<HlaAllele,int[]> alleleEntry : sourceEntry.getValue().entrySet())
                {
                    final HlaAllele allele = alleleEntry.getKey();
                    final int[] counts = alleleEntry.getValue();

                    int total = counts[SHARED] + counts[SHARED_HLAY] + counts[UNIQUE];

                    writer.write(String.format("%s,%s,%d,%d,%d,%d",
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
            return;
        }
    }

    private void updateMiscCounts(final int[] miscCounts, final List<Integer> fragAminoAcidLoci)
    {
        if(fragAminoAcidLoci.contains(227))
            ++miscCounts[Y0101_X];

        int exon3Start = A_EXON_BOUNDARIES.get(2) + 1;
        int exon3End = A_EXON_BOUNDARIES.get(3);

        if(fragAminoAcidLoci.get(0) >= exon3Start && fragAminoAcidLoci.get(fragAminoAcidLoci.size() - 1) <= exon3End)
            ++miscCounts[EXON_3];
    }

    private void writeFragmentInfo(
            final Fragment fragment, final List<Integer> fragAminoAcidLoci, final List<HlaAllele> hlayAlleles, boolean isShared)
    {
        try
        {
            if(mWriter == null)
            {
                String fileName = mOutputPrefix + ".hlay.fragments.csv";
                mWriter = createBufferedWriter(fileName, false);

                mWriter.write("ReadId,ReadInfo,Alleles,Shared,Loci");
                mWriter.newLine();
            }

            StringJoiner sjAlleles = new StringJoiner(ITEM_DELIM);
            hlayAlleles.forEach(x -> sjAlleles.add(x.toString()));

            StringJoiner sjLoci = new StringJoiner(ITEM_DELIM);
            fragAminoAcidLoci.forEach(x -> sjLoci.add(x.toString()));

            mWriter.write(String.format("%s,%s,%s,%s,%s",
                    fragment.id(), fragment.readInfo(), sjAlleles.toString(), isShared, sjLoci.toString()));
            mWriter.newLine();

        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to write HLA-Y fragments file: {}", e.toString());
            return;
        }
    }

}
