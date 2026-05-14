package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

// per-chromosome interval set covering ribosomal RNA contamination zones. Combines three sources:
//   (a) ensembl genes whose transcripts carry BioType in {rRNA, Mt_rRNA}, gene range [GeneStart, GeneEnd]
//   (b) hardcoded V38 acrocentric short-arm ranges (chr13/14/15/21/22 p-arms) which carry the rDNA
//       repeat arrays. These are biological constants and aren't fully captured by ensembl annotation.
//   (c) wholesale exclusion of contigs whose name starts with chr22_ or chrUn_ (V38). The chr22 random
//       contig and the chrUn decoys are dominated by rDNA junk — observed empirically as the top sinks
//       for rRNA reads. Any read mapping to a contig with these prefixes is treated as rRNA.
// (b) and (c) are V38-only; on V37 the index falls back to biotype-only with a warn.
// Lookups are inclusive at both ends; overlapping intervals are merged.
public class RrnaRegionIndex
{
    public static final String BIOTYPE_RRNA = "rRNA";
    public static final String BIOTYPE_MT_RRNA = "Mt_rRNA";

    // V38 acrocentric p-arm ranges. The short arms of these chromosomes are nucleolar organizer
    // regions carrying the 45S rDNA tandem arrays — STAR / bwa dump ambiguous rRNA reads here even
    // when ensembl annotation doesn't cover the locus. End coords sit just below the centromere start
    // so we capture the full rDNA-bearing p-arm without spilling into the centromere proper.
    private static final Interval[] V38_ACROCENTRIC_P_ARMS = {
            new Interval(1, 16_000_000),  // chr13
            new Interval(1, 16_000_000),  // chr14
            new Interval(1, 17_000_000),  // chr15
            new Interval(1, 12_000_000),  // chr21
            new Interval(1, 13_500_000),  // chr22
    };

    private static final String[] V38_ACROCENTRIC_CHROMOSOMES = {"chr13", "chr14", "chr15", "chr21", "chr22"};

    // contig-name prefixes for V38 unplaced / random contigs that are dominated by rDNA junk. Any read
    // mapping to a contig matching one of these is treated as rRNA regardless of position.
    private static final String[] V38_EXCLUDED_CONTIG_PREFIXES = {"chr22_", "chrUn_"};

    private final Map<String, List<Interval>> mIntervalsByChromosome;
    private final String[] mExcludedContigPrefixes;
    private final int mGeneCount;

    private RrnaRegionIndex(
            final Map<String, List<Interval>> intervalsByChromosome, final String[] excludedContigPrefixes, final int geneCount)
    {
        mIntervalsByChromosome = intervalsByChromosome;
        mExcludedContigPrefixes = excludedContigPrefixes;
        mGeneCount = geneCount;
    }

    public int geneCount()
    {
        return mGeneCount;
    }

    // returns true if the contig is wholly excluded by prefix, or if [start, end] (inclusive) overlaps
    // any rRNA region on the given chromosome.
    public boolean overlaps(final String chromosome, final int start, final int end)
    {
        for(final String prefix : mExcludedContigPrefixes)
        {
            if(chromosome.startsWith(prefix))
                return true;
        }

        final List<Interval> intervals = mIntervalsByChromosome.get(chromosome);
        if(intervals == null || intervals.isEmpty())
            return false;

        // binary search for the rightmost interval with Start <= end
        int lo = 0;
        int hi = intervals.size() - 1;
        int candidate = -1;
        while(lo <= hi)
        {
            int mid = (lo + hi) >>> 1;
            if(intervals.get(mid).Start <= end)
            {
                candidate = mid;
                lo = mid + 1;
            }
            else
            {
                hi = mid - 1;
            }
        }

        if(candidate < 0)
            return false;

        // walk back over any intervals whose End reaches our start. After merge there's typically
        // a single overlapping interval, but unmerged overlapping rRNA genes are tolerated.
        for(int i = candidate; i >= 0; --i)
        {
            final Interval interval = intervals.get(i);
            if(interval.End < start)
                return false;
            if(interval.Start <= end && interval.End >= start)
                return true;
        }
        return false;
    }

    public static RrnaRegionIndex build(final ConfigBuilder configBuilder, final RefGenomeVersion refGenomeVersion)
    {
        final EnsemblDataCache cache = new EnsemblDataCache(configBuilder);
        cache.setRequiredData(false, false, false, false);
        cache.load(false);
        return build(cache, refGenomeVersion);
    }

    static RrnaRegionIndex build(final EnsemblDataCache cache, final RefGenomeVersion refGenomeVersion)
    {
        final Map<String, List<TranscriptData>> transcriptsByGene = cache.getTranscriptDataMap();
        final Map<String, List<GeneData>> genesByChromosome = cache.getChrGeneDataMap();

        final Set<String> rrnaGeneIds = new HashSet<>();
        for(final Map.Entry<String, List<TranscriptData>> entry : transcriptsByGene.entrySet())
        {
            for(final TranscriptData transcript : entry.getValue())
            {
                if(isRrnaBiotype(transcript.BioType))
                {
                    rrnaGeneIds.add(entry.getKey());
                    break;
                }
            }
        }

        final Map<String, List<Interval>> rawIntervalsByChromosome = new HashMap<>();
        for(final Map.Entry<String, List<GeneData>> entry : genesByChromosome.entrySet())
        {
            final List<Interval> intervals = new ArrayList<>();
            for(final GeneData gene : entry.getValue())
            {
                if(rrnaGeneIds.contains(gene.GeneId))
                    intervals.add(new Interval(gene.GeneStart, gene.GeneEnd));
            }
            if(!intervals.isEmpty())
                rawIntervalsByChromosome.put(entry.getKey(), intervals);
        }

        final String[] excludedContigPrefixes;
        if(refGenomeVersion == RefGenomeVersion.V38)
        {
            for(int i = 0; i < V38_ACROCENTRIC_CHROMOSOMES.length; ++i)
            {
                rawIntervalsByChromosome
                        .computeIfAbsent(V38_ACROCENTRIC_CHROMOSOMES[i], k -> new ArrayList<>())
                        .add(V38_ACROCENTRIC_P_ARMS[i]);
            }
            excludedContigPrefixes = V38_EXCLUDED_CONTIG_PREFIXES;
        }
        else
        {
            RD_LOGGER.warn("rRNA filter: acrocentric p-arm and unplaced-contig exclusions are V38-only; "
                    + "running on ref genome version {} with biotype-only coverage", refGenomeVersion);
            excludedContigPrefixes = new String[0];
        }

        final Map<String, List<Interval>> intervalsByChromosome = new HashMap<>();
        for(final Map.Entry<String, List<Interval>> entry : rawIntervalsByChromosome.entrySet())
        {
            final List<Interval> intervals = entry.getValue();
            intervals.sort(Comparator.comparingInt(i -> i.Start));
            intervalsByChromosome.put(entry.getKey(), mergeOverlapping(intervals));
        }

        final int geneCount = rrnaGeneIds.size();
        final int chromosomeCount = intervalsByChromosome.size();
        RD_LOGGER.info("loaded {} rRNA gene region(s) across {} chromosome(s); excluded contig prefixes: {}",
                geneCount, chromosomeCount, excludedContigPrefixes.length == 0 ? "none" : String.join(",", excludedContigPrefixes));

        return new RrnaRegionIndex(intervalsByChromosome, excludedContigPrefixes, geneCount);
    }

    static boolean isRrnaBiotype(final String biotype)
    {
        return BIOTYPE_RRNA.equals(biotype) || BIOTYPE_MT_RRNA.equals(biotype);
    }

    private static List<Interval> mergeOverlapping(final List<Interval> sorted)
    {
        final List<Interval> merged = new ArrayList<>(sorted.size());
        Interval current = null;
        for(final Interval interval : sorted)
        {
            if(current == null)
            {
                current = interval;
                continue;
            }
            if(interval.Start <= current.End + 1)
            {
                current = new Interval(current.Start, Math.max(current.End, interval.End));
            }
            else
            {
                merged.add(current);
                current = interval;
            }
        }
        if(current != null)
            merged.add(current);
        return merged;
    }

    static final class Interval
    {
        final int Start;
        final int End;

        Interval(final int start, final int end)
        {
            Start = start;
            End = end;
        }
    }
}
