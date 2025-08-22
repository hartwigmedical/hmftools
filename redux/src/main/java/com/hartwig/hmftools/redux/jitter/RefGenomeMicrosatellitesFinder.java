package com.hartwig.hmftools.redux.jitter;

import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.GC_PROFILE;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.GC_PROFILE_DESC;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.loadChrGcProfileMap;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import static htsjdk.samtools.util.SequenceUtil.N;

import java.io.File;
import java.io.IOException;
import java.time.Duration;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.gc.ImmutableGCProfile;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.commons.lang3.Validate;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.StringUtil;

public class RefGenomeMicrosatellitesFinder
{
    public static final Logger MSI_LOGGER = LogManager.getLogger(RefGenomeMicrosatellitesFinder.class);

    private static final int BED_REGION_EXPANSION = 950;
    private static final int TARGET_SITE_COUNT = 3000;

    private static final double MIN_MAPPABILITY = 0.8;

    private static class Config
    {
        public final String RefGenomeFile;
        public final String OutputDir;

        public final RefGenomeVersion RefGenVersion;
        public final String GcProfilePath;
        public final String BedFile;

        public final int Threads;

        public Config(final ConfigBuilder configBuilder)
        {
            RefGenomeFile = configBuilder.getValue(REF_GENOME);
            OutputDir = parseOutputDir(configBuilder);
            RefGenVersion = RefGenomeVersion.from(configBuilder);
            GcProfilePath = configBuilder.getValue(GC_PROFILE);
            BedFile = configBuilder.getValue("bed_file");
            Threads = Math.max(parseThreads(configBuilder), 1);
        }

        public static void registerConfig(final ConfigBuilder configBuilder)
        {
            addRefGenomeVersion(configBuilder);
            configBuilder.addPath(REF_GENOME, true, REF_GENOME_CFG_DESC + ", required when using CRAM files");
            configBuilder.addPath(GC_PROFILE, true, GC_PROFILE_DESC);
            configBuilder.addPath("bed_file", false, "exta bed file input");

            addOutputDir(configBuilder);

            addLoggingOptions(configBuilder);
            addThreadOptions(configBuilder);

            addSpecificChromosomesRegionsConfig(configBuilder);
        }

        public boolean isValid()
        {
            checkCreateOutputDir(OutputDir);
            return true;
        }
    }

    static class Candidate
    {
        int startIndex;

        // where we are up to, it is also same as end
        int currentEndIndex;
        byte[] pattern;

        String patternString;

        boolean complete = false;

        Candidate(byte[] pattern, int startIndex, int endIndex)
        {
            this.pattern = pattern;
            this.startIndex = startIndex;
            this.currentEndIndex = endIndex;
            patternString = StringUtil.bytesToString(pattern);
        }

        byte nextBase() { return pattern[(currentEndIndex - startIndex) % pattern.length]; }

        int length() { return currentEndIndex - startIndex; }

        int numFullUnits() { return length() / pattern.length; }
    }

    // The jitter model merge some units together
    // Together 7 types
    enum UnitKey
    {
        _3_TO_5_BP("3-5pb"),
        A_T("A/T"),
        C_G("C/G"),
        AT_TA("AT/TA"),
        AG_GA_CT_TC("AG/GA/CT/TC"),
        AC_CA_GT_TG("AC/CA/GT/TG"),
        CG_GC("CG/GC");

        private final String unitKey;

        public String getUnitKey() { return unitKey; }

        UnitKey(String unitKey)
        {
            this.unitKey = unitKey;
        }

        public static UnitKey fromUnit(String unit)
        {
            if(unit.length() >= 3 && unit.length() <= 5)
            {
                return _3_TO_5_BP;
            }
            if(unit.equals("A") || unit.equals("T"))
            {
                return A_T;
            }
            if(unit.equals("C") || unit.equals("G"))
            {
                return C_G;
            }
            if(unit.equals("AT") || unit.equals("TA"))
            {
                return AT_TA;
            }
            if(unit.equals("AG") || unit.equals("GA") || unit.equals("CT") || unit.equals("TC"))
            {
                return AG_GA_CT_TC;
            }
            if(unit.equals("AC") || unit.equals("CA") || unit.equals("GT") || unit.equals("TG"))
            {
                return AC_CA_GT_TG;
            }
            if(unit.equals("CG") || unit.equals("GC"))
            {
                return CG_GC;
            }

            throw new IllegalStateException("unable to convert unit key: " + unit);
        }
    }

    private static class UnitRepeatKey
    {
        public final UnitKey Key;
        public final int NumRepeats;

        public UnitRepeatKey(final UnitKey unitKey, final int numRepeats)
        {
            Key = unitKey;
            NumRepeats = numRepeats;
        }

        @Override
        public boolean equals(final Object o)
        {
            if(this == o)
            {
                return true;
            }
            if(!(o instanceof UnitRepeatKey))
            {
                return false;
            }

            final UnitRepeatKey that = (UnitRepeatKey) o;

            if(NumRepeats != that.NumRepeats)
            {
                return false;
            }
            return Key == that.Key;
        }

        @Override
        public int hashCode()
        {
            int result = Key.hashCode();
            result = 31 * result + NumRepeats;
            return result;
        }

        @Override
        public String toString()
        {
            return Key.getUnitKey() + " x " + NumRepeats;
        }
    }

    private final Config mConfig;
    private final  Map<String,List<GCProfile>> mGcProfiles;
    private final Multimap<UnitRepeatKey, RefGenomeMicrosatellite> mAllMicrosatelliteSites = ArrayListMultimap.create();
    
    // this needs to be thread safe
    private final Multimap<UnitRepeatKey, RefGenomeMicrosatellite> mDownSampledMicrosatelliteSites = Multimaps.synchronizedListMultimap(ArrayListMultimap.create());

    public RefGenomeMicrosatellitesFinder(final ConfigBuilder configBuilder) throws IOException
    {
        mConfig = new Config(configBuilder);
        mGcProfiles = loadChrGcProfileMap(mConfig.GcProfilePath);
    }

    public int run() throws Exception
    {
        MSI_LOGGER.info("finding ms sites from ref genome: {}", mConfig.RefGenomeFile);
        Instant start = Instant.now();

        IndexedFastaSequenceFile refGenome = new IndexedFastaSequenceFile(new File(mConfig.RefGenomeFile));

        //
        RefGenomeMicrosatellitesFinder.findMicrosatellites(refGenome, JitterAnalyserConstants.MIN_MICROSAT_UNIT_COUNT,
                r -> {
                    populateMappability(r, mGcProfiles);

                    // put all into multimap
                    mAllMicrosatelliteSites.put(new UnitRepeatKey(UnitKey.fromUnit(r.unitString()), r.RepeatCount), r);
                });

        // now do downsampling
        downsampleSites(TARGET_SITE_COUNT);

        // now write only the downsampled sites to file
        String outputFile = RefGenomeMicrosatelliteFile.generateFilename(mConfig.OutputDir, mConfig.RefGenVersion);
        try (RefGenomeMicrosatelliteFile refGenomeMicrosatelliteFile = new RefGenomeMicrosatelliteFile(outputFile))
        {
            mDownSampledMicrosatelliteSites.values().forEach(refGenomeMicrosatelliteFile::writeRow);
        }

        MSI_LOGGER.info("wrote {} microsatellite sites into {}", mDownSampledMicrosatelliteSites.size(), outputFile);

        Duration duration = Duration.between(start, Instant.now());

        MSI_LOGGER.info("ref genome microsatellites finder took {}m {}s", duration.toMinutes(), duration.toSecondsPart());

        return 0;
    }

    // check if this is a valid repeat unit, i.e.
    // ATAT is not because it is actually 2 x AT
    static boolean isValidUnit(byte[] unit)
    {
        // check each subunit length to see if it is the same one
        for(int subunitLength = 1; subunitLength <= unit.length / 2; ++subunitLength)
        {
            boolean allMatch = true;

            if((unit.length % subunitLength) == 0)
            {
                for(int i = 0; i < subunitLength && allMatch; ++i)
                {
                    byte base = unit[i];
                    for(int j = subunitLength + i; j < unit.length; j += subunitLength)
                    {
                        if(unit[j] != base)
                        {
                            allMatch = false;
                            break;
                        }
                    }
                }
            }
            else
            {
                allMatch = false;
            }

            if (allMatch)
            {
                // found a subunit that matches the whole unit
                return false;
            }
        }
        // did not find a subunit
        return true;
    }

    public static void findMicrosatellites(
            final ReferenceSequenceFile referenceSequenceFile, int minRepeatLength,
            final Consumer<RefGenomeMicrosatellite> refGenomeMsConsumer)
    {
        int chunkSize = 100_000;
        findMicrosatellites(referenceSequenceFile, minRepeatLength, refGenomeMsConsumer, chunkSize);
    }

    // algorithm to find short tandem repeats
    // at each base, we try find longest candidate starting from this base.
    // if a microsatellite is found, we start again from the base after.
    //
    static void findMicrosatellites(
            final ReferenceSequenceFile referenceSequenceFile, int minNumRepeats,
            final Consumer<RefGenomeMicrosatellite> refGenomeMsConsumer, int chunkSize)
    {
        MutableInt microsatelliteCounter = new MutableInt(0);

        List<SAMSequenceRecord> seqRecords = referenceSequenceFile.getSequenceDictionary().getSequences()
                .stream()
                .filter(o -> HumanChromosome.contains(o.getContig()))
                .collect(Collectors.toList());

        for(SAMSequenceRecord sequenceRecord : seqRecords)
        {
            MSI_LOGGER.info("start processing chromosome {}", sequenceRecord.getContig());

            int length = sequenceRecord.getSequenceLength();

            // pending candidate, we do not accept a candidate when it is completed. We want to avoid
            // candidates that are too close to each other. They are dropped if too close.
            List<RefGenomeMicrosatellite> pendingMicrosatellies = new ArrayList<>();

            // current best candidate
            Candidate bestCandidate = null;

            // candidates
            List<Candidate> currentCandidates = new ArrayList<>();

            for(int start = 1; start <= length; start = start + chunkSize)
            {
                int endInclusive = Math.min(start + chunkSize, length) - 1;
                byte[] seq = referenceSequenceFile.getSubsequenceAt(sequenceRecord.getContig(), start, endInclusive).getBases();

                for(int i = 0; i < seq.length; ++i)
                {
                    byte base = seq[i];

                    if(base == N)
                    {
                        // if we hit an N, we clear everything
                        currentCandidates.clear();
                        bestCandidate = null;
                        continue;
                    }

                    // check all current candidates
                    ListIterator<Candidate> itr = currentCandidates.listIterator();
                    while(itr.hasNext())
                    {
                        Candidate c = itr.next();

                        // check if pattern is still valid
                        if(c.nextBase() == base)
                        {
                            // pattern continues
                            c.currentEndIndex++;
                            if(bestCandidate == null || bestCandidate.numFullUnits() < c.numFullUnits())
                            {
                                bestCandidate = c;
                            }
                        }
                        else
                        {
                            c.complete = true;
                            itr.remove();
                        }
                    }

                    // if the best candidate is completed, we want to check if it is a candidate
                    // we want to include
                    if(bestCandidate != null && bestCandidate.complete)
                    {
                        // check if the best candidate that is completed is a good one
                        int unitRepeatCount = bestCandidate.numFullUnits();

                        if(unitRepeatCount >= minNumRepeats)
                        {
                            // NOTE: for microsatellites with a pattern longer than one, we only accept full repeats
                            // i.e. ATGATGATGATGAT contains 4 full ATG plus AT at the end, even though AT is partial unit,
                            // they are excluded.

                            int baseLength = unitRepeatCount * bestCandidate.pattern.length;

                            // this is a microsatellite
                            RefGenomeMicrosatellite refGenomeMicrosatellite = new RefGenomeMicrosatellite(sequenceRecord.getContig(),
                                    bestCandidate.startIndex,
                                    bestCandidate.startIndex + baseLength - 1, // change to inclusive
                                    bestCandidate.pattern);
                            pendingMicrosatellies.add(refGenomeMicrosatellite);

                            // check the panding microsatellites, see if any can be accepted
                            checkPendingMicrosatellites(pendingMicrosatellies, refGenomeMsConsumer, microsatelliteCounter);
                        }

                        bestCandidate = null;
                    }

                    // also start new Candidates at this location
                    // all these candidates have the current base as the last base of the pattern
                    for(int j = Math.max(i - JitterAnalyserConstants.MAX_MICROSAT_UNIT_LENGTH + 1, 0); j <= i; ++j)
                    {
                        byte[] repeatUnit = Arrays.copyOfRange(seq, j, i + 1);

                        // check that this is a valid repeat unit, i.e. it is not multiple smaller unit
                        if(isValidUnit(repeatUnit))
                        {
                            // also check against existing patterns
                            if(currentCandidates.stream().noneMatch(o -> Arrays.equals(o.pattern, repeatUnit)))
                            {
                                Candidate c = new Candidate(repeatUnit, start + j, start + i + 1);
                                currentCandidates.add(c);
                            }
                        }
                    }
                }
            }

            if(pendingMicrosatellies.size() == 1)
            {
                microsatelliteCounter.increment();
                refGenomeMsConsumer.accept(pendingMicrosatellies.get(0));
            }

            MSI_LOGGER.info("finished chromosome {}", sequenceRecord.getSequenceName());
        }

        MSI_LOGGER.info("found {} microsatellite regions in ref genome", microsatelliteCounter);
    }

    // the aim of this code is to remove microsatellites that are too close to each other
    // We do not try to merge them for now, even though tools such as MsDetector would.
    // i.e. AAAAATAAAAAAA
    static void checkPendingMicrosatellites(
            final List<RefGenomeMicrosatellite> pendingMicrosatellies, final Consumer<RefGenomeMicrosatellite> refGenomeMsConsumer,
            final MutableInt microsatelliteCounter)
    {
        int groupStart = 0;

        // check the panding microsatellites, see if any can be accepted
        for(int i = 0; i < pendingMicrosatellies.size() - 1; ++i)
        {
            RefGenomeMicrosatellite ms1 = pendingMicrosatellies.get(i);
            RefGenomeMicrosatellite ms2 = pendingMicrosatellies.get(i + 1);

            if((ms2.referenceStart() - ms1.referenceEnd()) > JitterAnalyserConstants.MIN_ADJACENT_MICROSAT_DISTANCE)
            {
                // previous group finished, if previous group only has 1 ms, we accept it, otherwise
                // remove them all from the pending list
                if(i == groupStart)
                {
                    // only 1 item, accept this
                    refGenomeMsConsumer.accept(ms1);
                    microsatelliteCounter.increment();
                    MSI_LOGGER.trace("microsatellite: {}", ms1);
                }
                else
                {
                    // ms are too close to each other, remove them
                    for(int j = groupStart; j <= i; ++j)
                    {
                        MSI_LOGGER.trace("reject microsatellite as too close to neighbour: {}", pendingMicrosatellies.get(j));
                    }
                }

                // update group start
                groupStart = i + 1;
            }
        }

        // we can remove anything before groupStart
        if(groupStart > 0)
        {
            pendingMicrosatellies.subList(0, groupStart).clear();
        }
    }

    private static void populateMappability(
            final RefGenomeMicrosatellite refGenomeMicrosatellite, final Map<String,List<GCProfile>> gcProfileMap)
    {
        List<GCProfile> gcProfiles = gcProfileMap.get(refGenomeMicrosatellite.chromosome());
        int siteMid = refGenomeMicrosatellite.Region.start() + (refGenomeMicrosatellite.Region.baseLength()) / 2;

        GCProfile endKey = ImmutableGCProfile.builder().from(gcProfiles.get(0)).end(siteMid).build();

        // find the position using binary search
        int index = Collections.binarySearch(gcProfiles, endKey, Comparator.comparingInt(GCProfile::end));

        // we go backwards and forwards
        if(index < 0)
        {
            index = -index - 1;
        }

        if(index < gcProfiles.size())
        {
            GCProfile gcProfile = gcProfiles.get(index);
            Validate.isTrue(gcProfile.overlaps(refGenomeMicrosatellite.Region.genomeRegion()));
            refGenomeMicrosatellite.setMappability(gcProfile.mappablePercentage());
        }
        else
        {
            refGenomeMicrosatellite.setMappability(0);
            MSI_LOGGER.warn("microsatellite site({}) gc profile not found", refGenomeMicrosatellite.Region);
        }
    }

    // filter the microsatellites such that each type of (unit, length) is approximately the target count
    void downsampleSites(int targetCountPerType) throws ExecutionException, InterruptedException
    {
        MSI_LOGGER.info("filtering microsatellite sites, target count per type = {}", targetCountPerType);

        Map<String,List<BaseRegion>> regionsToKeep = ChrBaseRegion.loadChrBaseRegions(mConfig.BedFile);


        int numDigits = Integer.toString(mConfig.Threads - 1).length();
        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("thread-%0" + numDigits + "d").build();
        ExecutorService executorService = Executors.newFixedThreadPool(mConfig.Threads, namedThreadFactory);

        final List<CompletableFuture<Void>> futures = new ArrayList<>();

        // first we need to decide how much to downsample by. We do so by working out how many sites need to be
        for(UnitRepeatKey unitRepeatKey : mAllMicrosatelliteSites.keySet())
        {
            Runnable task = () -> downSampleMsType(targetCountPerType, regionsToKeep, unitRepeatKey);
            futures.add(CompletableFuture.runAsync(task, executorService));
        }

        // wait for completion
        CompletableFuture.allOf(futures.toArray(CompletableFuture[]::new)).get();

        MSI_LOGGER.info("filtered {} microsatellite sites down to {}",
                mAllMicrosatelliteSites.size(), mDownSampledMicrosatelliteSites.size());
    }

    private void downSampleMsType(
            final int targetCountPerType, final Map<String,List<BaseRegion>> regionsToKeep,
            final UnitRepeatKey unitRepeatKey)
    {
        Collection<RefGenomeMicrosatellite> filteredList = mDownSampledMicrosatelliteSites.get(unitRepeatKey);
        Collection<RefGenomeMicrosatellite> allList = mAllMicrosatelliteSites.get(unitRepeatKey);

        MSI_LOGGER.info("[{}] filtering {} microsatellite sites", unitRepeatKey, allList.size());

        if(allList.size() <= targetCountPerType)
        {
            // no filtering for this ms type
            filteredList.addAll(allList);
            return;
        }

        List<RefGenomeMicrosatellite> sitesOutsideBedRegions = new ArrayList<>();

        // first work out which are retained by the bed regions
        // we do it first by taking all inside the bed regions + 1000 bases around those bed regions
        // if there are too many, then we only take those inside the bed regions.

        // As a speed optimisation, we assume that if we have more than 100x more sites than what we need
        // then we should use the smaller bed regions and skip the expanded bed
        int bedRegionExpansion = allList.size() > 100 * targetCountPerType ? 0 : BED_REGION_EXPANSION;
        double mappabilityCutoff = MIN_MAPPABILITY;

        while(true)
        {
            for(RefGenomeMicrosatellite r : allList)
            {
                if(r.mappability() < mappabilityCutoff)
                {
                    continue;
                }

                boolean isInBed = false;
                for(BaseRegion region : regionsToKeep.get(r.chromosome()))
                {
                    if(positionsOverlap(
                            r.Region.start(), r.Region.end(),
                            region.start() - bedRegionExpansion, region.end() + bedRegionExpansion))
                    {
                        isInBed = true;
                        break;
                    }
                }
                if(isInBed)
                {
                    filteredList.add(r);
                }
                else
                {
                    sitesOutsideBedRegions.add(r);
                }
            }

            if(bedRegionExpansion > 0 && filteredList.size() > targetCountPerType * 2)
            {
                // if we got way too many sites inside bed region try again with no expanded region and also mappability cutoff of 1.0
                MSI_LOGGER.debug("[{}] too many sites in bed + expanded region: {}, try again with no expanded region",
                        unitRepeatKey, filteredList.size());
                bedRegionExpansion = 0;
                mappabilityCutoff = 1.0;
                filteredList.clear();
                sitesOutsideBedRegions.clear();
                continue;
            }
            break;
        }

        // now decide how many to filter out from the rest
        int numSitesInBed = filteredList.size();

        double frac = ((double) targetCountPerType - numSitesInBed) / sitesOutsideBedRegions.size();
        if(frac > 0.0)
        {
            if(frac < 1.0)
            {
                final Random random = new Random();
                sitesOutsideBedRegions.stream().filter(o -> random.nextDouble() <= frac).forEach(filteredList::add);
            }
            else
            {
                filteredList.addAll(sitesOutsideBedRegions);
            }
        }

        MSI_LOGGER.info("[{}] filtered {} microsatellite sites down to {}, sites in bed: {}, sites outside bed: {}",
                unitRepeatKey, allList.size(), filteredList.size(), numSitesInBed, sitesOutsideBedRegions.size());
    }

    public static void main(final String... args) throws Exception
    {
        ConfigBuilder configBuilder = new ConfigBuilder("RefGenomeMicrosatellitesFinder");
        Config.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        System.exit(new RefGenomeMicrosatellitesFinder(configBuilder).run());
    }
}
