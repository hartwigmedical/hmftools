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
import java.util.NavigableSet;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multiset;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.utils.AminoAcid;
import com.hartwig.hmftools.lilac.utils.Nucleotide;

import org.jetbrains.annotations.Nullable;

public final class SequenceCount
{
    private final int mMinFilterDepth;
    private final double mMinEvidenceFactor;

    private final NavigableMap<Integer, Multiset<String>> mSeqCountsByLoci;

    private SequenceCount(final int minFilterDepth, final double minEvidenceFactor, final NavigableMap<Integer, Multiset<String>> seqCounts)
    {
        mMinFilterDepth = minFilterDepth;
        mMinEvidenceFactor = minEvidenceFactor;
        mSeqCountsByLoci = seqCounts;
    }

    @VisibleForTesting
    public SequenceCount(final Map<String, Integer>[] seqCounts)
    {
        mSeqCountsByLoci = Maps.newTreeMap();
        for(int locus = 0; locus < seqCounts.length; locus++)
        {
            Multiset<String> seqCount = HashMultiset.create();
            for(Map.Entry<String, Integer> entry : seqCounts[locus].entrySet())
            {
                seqCount.add(entry.getKey(), entry.getValue());
            }

            mSeqCountsByLoci.put(locus, seqCount);
        }

        mMinFilterDepth = Integer.MAX_VALUE;
        mMinEvidenceFactor = 0.0;
    }

    public static SequenceCount nucleotides(final int minFilterDepth, final double minEvidenceFactor, final Iterable<Fragment> fragments)
    {
        NavigableMap<Integer, Multiset<String>> seqCountsByLoci = Maps.newTreeMap();
        for(Fragment fragment : fragments)
        {
            for(Nucleotide nucleotide : fragment.nucleotidesByLoci().values())
            {
                int locus = nucleotide.locus();
                String bases = nucleotide.bases();
                seqCountsByLoci.computeIfAbsent(locus, k -> HashMultiset.create());
                seqCountsByLoci.get(locus).add(bases);
            }
        }

        return new SequenceCount(minFilterDepth, minEvidenceFactor, seqCountsByLoci);
    }

    public static SequenceCount aminoAcids(final int minFilterDepth, final double minEvidenceFactor, final Iterable<Fragment> fragments)
    {
        NavigableMap<Integer, Multiset<String>> seqCountsByLoci = Maps.newTreeMap();
        for(Fragment fragment : fragments)
        {
            for(AminoAcid aminoAcid : fragment.aminoAcidsByLoci().values())
            {
                int locus = aminoAcid.locus();
                String acid = aminoAcid.acid();
                seqCountsByLoci.computeIfAbsent(locus, k -> HashMultiset.create());
                seqCountsByLoci.get(locus).add(acid);
            }
        }

        return new SequenceCount(minFilterDepth, minEvidenceFactor, seqCountsByLoci);
    }

    public Multiset<String> get(final int locus)
    {
        return mSeqCountsByLoci.getOrDefault(locus, HashMultiset.create());
    }

    public NavigableMap<Integer, Multiset<String>> seqCountsByLoci()
    {
        return mSeqCountsByLoci;
    }

    public NavigableSet<Integer> heterozygousLoci()
    {
        return mSeqCountsByLoci.keySet().stream()
                .filter(this::isHeterozygous)
                .collect(Collectors.toCollection(Sets::newTreeSet));
    }

    public NavigableSet<Integer> homozygousLoci()
    {
        return mSeqCountsByLoci.keySet().stream()
                .filter(this::isHomozygous)
                .collect(Collectors.toCollection(Sets::newTreeSet));
    }

    private boolean isHomozygous(final int locus)
    {
        return getMinEvidenceSequences(locus).size() == 1;
    }

    private boolean isHeterozygous(final int locus)
    {
        return getMinEvidenceSequences(locus).size() > 1;
    }

    public List<String> getMinEvidenceSequences(final int locus)
    {
        return getMinEvidenceSequences(locus, null);
    }

    public List<String> getMinEvidenceSequences(final int locus, @Nullable final Double minEvidenceFactor)
    {
        Multiset<String> seqCounts = mSeqCountsByLoci.get(locus);
        if(seqCounts == null)
        {
            return Lists.newArrayList();
        }

        int support = seqCounts.size();
        if(support < mMinFilterDepth)
        {
            return seqCounts.entrySet().stream()
                    .map(Multiset.Entry::getElement)
                    .collect(Collectors.toList());
        }

        double factor = minEvidenceFactor == null ? mMinEvidenceFactor : minEvidenceFactor;
        int evidenceMinCount = (int) ceil(support * factor);
        return seqCounts.entrySet().stream()
                .filter(x -> x.getCount() >= evidenceMinCount)
                .map(Multiset.Entry::getElement)
                .collect(Collectors.toList());
    }

    public String getMaxCountSequence(final int locus)
    {
        Multiset<String> seqCounts = mSeqCountsByLoci.get(locus);
        if(seqCounts == null || seqCounts.isEmpty())
        {
            return "-";
        }

        Multiset.Entry<String> maxEntry = seqCounts.entrySet().stream().max(Comparator.comparingInt(Multiset.Entry::getCount)).orElse(null);
        return maxEntry.getElement();
    }

    public int depth(final int locus)
    {
        return mSeqCountsByLoci.get(locus).size();
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
                {
                    continue;
                }

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
                {
                    continue;
                }

                Set<String> combinedSeqs = Sets.newHashSet();
                combinedSeqs.addAll(aSeqs);
                combinedSeqs.addAll(bSeqs);
                aHetLociMap.put(locus, combinedSeqs);
                bHetLociMap.put(locus, combinedSeqs);
            }
        }

        return geneHetLociMap;
    }

    private Map<Integer, Set<String>> extractHeterozygousLociSequences(final Iterable<HlaSequenceLoci> extraSequences)
    {
        NavigableMap<Integer, Set<String>> lociSeqMap = Maps.newTreeMap();
        for(Map.Entry<Integer, Multiset<String>> entry : mSeqCountsByLoci.entrySet())
        {
            int locus = entry.getKey();
            Set<String> aminoAcids = Sets.newHashSet(getMinEvidenceSequences(locus));
            for(HlaSequenceLoci extraSeqLoci : extraSequences)
            {
                if(locus >= extraSeqLoci.length())
                {
                    continue;
                }

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

            for(Map.Entry<Integer, Multiset<String>> seqCountsEntry : mSeqCountsByLoci.entrySet())
            {
                int locus = seqCountsEntry.getKey();
                Multiset<String> seqCounts = seqCountsEntry.getValue();

                StringJoiner lineBuilder = new StringJoiner("\t");
                lineBuilder.add(String.valueOf(locus));

                Iterator<Multiset.Entry<String>> sortedCounts = seqCounts.entrySet().stream()
                        .sorted(Comparator.comparingInt((Multiset.Entry<String> x) -> x.getCount()).reversed())
                        .limit(5)
                        .iterator();

                while(sortedCounts.hasNext())
                {
                    Multiset.Entry<String> count = sortedCounts.next();
                    lineBuilder.add(count.getElement());
                    lineBuilder.add(String.valueOf(count.getCount()));
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
