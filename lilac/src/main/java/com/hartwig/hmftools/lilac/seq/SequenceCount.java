package com.hartwig.hmftools.lilac.seq;

import static java.lang.Math.ceil;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.lilac.GeneCache.shortGeneName;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_A;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_B;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_C;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.utils.AminoAcid;
import com.hartwig.hmftools.lilac.utils.Nucleotide;

import org.jetbrains.annotations.Nullable;

public final class SequenceCount
{
    private final int mMinCount;

    private final NavigableMap<Integer, Map<String, Integer>> mSeqCountsByLoci;

    @VisibleForTesting
    public SequenceCount(int minCount, final NavigableMap<Integer, Map<String, Integer>> seqCounts)
    {
        mMinCount = minCount;
        mSeqCountsByLoci = seqCounts;
    }

    @VisibleForTesting
    public SequenceCount(int minCount, final Map<String, Integer>[] seqCounts)
    {
        mMinCount = minCount;
        mSeqCountsByLoci = Maps.newTreeMap();
        for(int i = 0; i < seqCounts.length; i++)
            mSeqCountsByLoci.put(i, seqCounts[i]);
    }

    public Map<String, Integer> get(int locus)
    {
        Map<String, Integer> seqCounts = mSeqCountsByLoci.get(locus);
        if(seqCounts == null)
        {
            LL_LOGGER.error("invalid sequence count index({}) look-up: loci({})", locus, Lists.newArrayList(mSeqCountsByLoci.keySet()));
            return Maps.newHashMap();
        }

        return seqCounts;
    }

    public static int calcMinCount(int minFilterDepth, double minFactor, int depth)
    {
        return depth < minFilterDepth ? 0 : (int) ceil(minFactor * depth);
    }

    public static SequenceCount nucleotides(int minFilterDepth, double minFactor, final Collection<Fragment> fragments)
    {
        NavigableMap<Integer, Map<String, Integer>> seqCountsByLoci = Maps.newTreeMap();
        for(Fragment fragment : fragments)
        {
            for(Nucleotide nucleotide : fragment.nucleotidesByLoci().values())
            {
                int locus = nucleotide.locus();
                String bases = nucleotide.bases();
                increment(seqCountsByLoci, locus, bases);
            }
        }

        return new SequenceCount(calcMinCount(minFilterDepth, minFactor, fragments.size()), seqCountsByLoci);
    }

    public static SequenceCount aminoAcids(int minFilterDepth, double minFactor, final Collection<Fragment> fragments)
    {

        NavigableMap<Integer, Map<String, Integer>> seqCountsByLoci = Maps.newTreeMap();
        for(Fragment fragment : fragments)
        {
            for(AminoAcid aminoAcid : fragment.aminoAcidsByLoci().values())
            {
                int locus = aminoAcid.locus();
                String acid = aminoAcid.acid();
                increment(seqCountsByLoci, locus, acid);
            }
        }

        return new SequenceCount(calcMinCount(minFilterDepth, minFactor, fragments.size()), seqCountsByLoci);
    }

    public NavigableMap<Integer, Map<String, Integer>> seqCountsByLoci() { return mSeqCountsByLoci; }

    public List<Integer> heterozygousLoci()
    {
        List<Integer> loci = Lists.newArrayList();
        for(int locus : mSeqCountsByLoci.keySet())
        {
            if(isHeterozygous(locus))
                loci.add(locus);
        }

        return loci;
    }

    public List<Integer> homozygousLoci()
    {
        List<Integer> loci = Lists.newArrayList();
        for(int locus : mSeqCountsByLoci.keySet())
        {
            if(isHomozygous(locus))
                loci.add(locus);
        }

        return loci;
    }

    private boolean isHomozygous(int locus)
    {
        Map<String, Integer> seqCounts = get(locus);
        return seqCounts.values().stream().filter(x -> x >= mMinCount).count() == 1;
    }

    private boolean isHeterozygous(int locus)
    {
        Map<String, Integer> seqCounts = get(locus);
        return seqCounts.values().stream().filter(x -> x >= mMinCount).count() > 1;
    }

    public List<String> getMinCountOrVafSequences(int locus, @Nullable final Double minLocalVAF)
    {
        Map<String, Integer> seqCounts = mSeqCountsByLoci.get(locus);
        if(seqCounts == null)
            return Lists.newArrayList();

        double vafMinCount = 0;
        if(minLocalVAF != null)
        {
            int coverageAtCurrentPosition = seqCounts.values().stream().mapToInt(Integer::intValue).sum();
            vafMinCount = coverageAtCurrentPosition * minLocalVAF;
        }

        double minCount = minLocalVAF != null ? vafMinCount : mMinCount;
        return seqCounts.entrySet().stream()
                .filter(x -> x.getValue() >= minCount)
                .map(Map.Entry::getKey).collect(Collectors.toList());
    }

    public List<String> getMinCountSequences(int locus)
    {
        return getMinCountOrVafSequences(locus, null);
    }

    public String getMaxCountSequence(int locus)
    {
        Map<String, Integer> seqCounts = get(locus);

        String sequence = "-";
        int maxCountSeq = 0;

        for(Map.Entry<String,Integer> entry : seqCounts.entrySet())
        {
            if(entry.getValue() > maxCountSeq)
            {
                sequence = entry.getKey();
                maxCountSeq = entry.getValue();
            }
        }

        return sequence;
    }

    public int depth(int locus)
    {
        Map<String, Integer> seqCounts = get(locus);
        return seqCounts.values().stream().mapToInt(x -> x).sum();
    }

    private static void increment(final NavigableMap<Integer, Map<String, Integer>> seqCountsByLoci, int locus, final String aminoAcid)
    {
        seqCountsByLoci.computeIfAbsent(locus, l -> Maps.newTreeMap());
        Map<String, Integer> seqCounts = seqCountsByLoci.get(locus);
        seqCounts.merge(aminoAcid, 1, Integer::sum);
    }

    public static Map<String, Map<Integer, Set<String>>> extractHeterozygousLociSequences(final Map<String, SequenceCount> geneCountsMap,
            final Collection<HlaSequenceLoci> extraSeqLoci)
    {
        Map<String, Map<Integer, Set<String>>> geneHetLociMap = Maps.newHashMap();
        for(Map.Entry<String, SequenceCount> geneEntry : geneCountsMap.entrySet())
        {
            String gene = shortGeneName(geneEntry.getKey());
            SequenceCount sequenceCounts = geneEntry.getValue();
            List<HlaSequenceLoci> geneExtraSeqLoci = extraSeqLoci.stream().filter(x -> x.Allele.Gene.equals(gene)).toList();
            Map<Integer, Set<String>> hetLociMap = sequenceCounts.extractHeterozygousLociSequences(geneExtraSeqLoci);
            geneHetLociMap.put(gene, hetLociMap);
        }

        // for recovered alleles (the extra-seq-loci), any additional amino acid location prior to 337 needs to be evaluated against
        // all 3 genes and added to all of them. From 338 onwards, A and B should be shared with each other, but C needs to be separate.
        Map<Integer, Set<String>> aHetLociMap = geneHetLociMap.get(GENE_A);
        Map<Integer, Set<String>> bHetLociMap = geneHetLociMap.get(GENE_B);
        Map<Integer, Set<String>> cHetLociMap = geneHetLociMap.get(GENE_C);

        for(int locus = 0; locus <= 348; ++locus)
        {
            Set<String> aSeqs = aHetLociMap.containsKey(locus) ? aHetLociMap.get(locus) : Sets.newHashSet();
            Set<String> bSeqs = bHetLociMap.containsKey(locus) ? bHetLociMap.get(locus) : Sets.newHashSet();
            Set<String> cSeqs = cHetLociMap.containsKey(locus) ? cHetLociMap.get(locus) : Sets.newHashSet();

            if(locus <= 337)
            {
                if(aSeqs.isEmpty() && bSeqs.isEmpty() && cSeqs.isEmpty())
                    continue;

                Set<String> combinedSeqs = Sets.newHashSet();
                combinedSeqs.addAll(aSeqs);
                combinedSeqs.addAll(bSeqs);
                combinedSeqs.addAll(cSeqs);
                aHetLociMap.put(locus, combinedSeqs);
                bHetLociMap.put(locus, combinedSeqs);
                cHetLociMap.put(locus, combinedSeqs);
            }
            else
            {
                // only A and B
                if(aSeqs.isEmpty() && bSeqs.isEmpty())
                    continue;

                Set<String> combinedSeqs = Sets.newHashSet();
                combinedSeqs.addAll(aSeqs);
                combinedSeqs.addAll(bSeqs);
                aHetLociMap.put(locus, combinedSeqs);
                bHetLociMap.put(locus, combinedSeqs);
            }
        }

        return geneHetLociMap;
    }

    private Map<Integer, Set<String>> extractHeterozygousLociSequences(final List<HlaSequenceLoci> extraSequences)
    {
        NavigableMap<Integer, Set<String>> lociSeqMap = Maps.newTreeMap();
        for(Map.Entry<Integer, Map<String, Integer>> entry : mSeqCountsByLoci.entrySet())
        {
            int locus = entry.getKey();
            Map<String, Integer> seqCounts = entry.getValue();

            Set<String> aminoAcids = seqCounts.entrySet().stream()
                    .filter(x -> x.getValue() >= mMinCount).map(Map.Entry::getKey).collect(Collectors.toSet());

            for(HlaSequenceLoci extraSeqLoci : extraSequences)
            {
                if(locus >= extraSeqLoci.length())
                    continue;

                String lociSeq = extraSeqLoci.sequence(locus);
                aminoAcids.add(lociSeq);
            }

            if(aminoAcids.size() > 1)
            {
                lociSeqMap.put(locus, aminoAcids);
            }
        }

        return lociSeqMap;
    }

    public void writeVertically(final String fileName)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(fileName, false);

            for(Map.Entry<Integer, Map<String, Integer>> seqCountsEntry : mSeqCountsByLoci.entrySet())
            {
                int locus = seqCountsEntry.getKey();
                Map<String, Integer> seqMap = seqCountsEntry.getValue();

                StringJoiner lineBuilder = new StringJoiner("\t");
                lineBuilder.add(String.valueOf(locus));

                Iterator<Map.Entry<String, Integer>> sortedCounts = seqMap.entrySet().stream()
                        .sorted(Comparator.comparingInt((Map.Entry<String, Integer> x) -> x.getValue()).reversed())
                        .limit(5)
                        .iterator();

                while(sortedCounts.hasNext())
                {
                    Map.Entry<String, Integer> count = sortedCounts.next();
                    lineBuilder.add(count.getKey());
                    lineBuilder.add(String.valueOf(count.getValue()));
                }

                writer.write(lineBuilder.toString());
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to write {}: {}", fileName, e.toString());
        }
    }
}
