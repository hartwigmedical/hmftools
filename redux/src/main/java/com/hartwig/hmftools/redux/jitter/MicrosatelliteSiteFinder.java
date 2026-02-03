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
import static com.hartwig.hmftools.redux.ReduxConfig.APP_NAME;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

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

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.StringUtil;

public class MicrosatelliteSiteFinder
{
    private static final int BED_REGION_EXPANSION = 950;
    private static final int MIN_TARGET_SITE_COUNT = 15_000;

    private static final double MIN_MAPPABILITY = 0.7;

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
        public final int length;

        public String getUnitKey() { return unitKey; }

        UnitKey(String unitKey)
        {
            this.unitKey = unitKey;
            if(unitKey.contains("/"))
                this.length = unitKey.split("/")[0].length();
            else
                this.length = 3;
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
    private final Multimap<UnitRepeatKey, MicrosatelliteSite> mAllMicrosatelliteSites = ArrayListMultimap.create();
    
    // this needs to be thread safe
    private final Multimap<UnitRepeatKey, MicrosatelliteSite> mDownSampledMicrosatelliteSites = Multimaps.synchronizedListMultimap(ArrayListMultimap.create());

    public MicrosatelliteSiteFinder(final ConfigBuilder configBuilder) throws IOException
    {
        mConfig = new Config(configBuilder);
        mGcProfiles = loadChrGcProfileMap(mConfig.GcProfilePath);
    }

    public int run() throws Exception
    {
        RD_LOGGER.info("finding ms sites from ref genome: {}", mConfig.RefGenomeFile);
        Instant start = Instant.now();

        IndexedFastaSequenceFile refGenome = new IndexedFastaSequenceFile(new File(mConfig.RefGenomeFile));

        MicrosatelliteSiteFinder.findMicrosatellites(refGenome, JitterConstants.MIN_MICROSAT_UNIT_COUNT,
                r -> {
                    populateMappability(r, mGcProfiles);

                    // put all into multimap
                    mAllMicrosatelliteSites.put(new UnitRepeatKey(UnitKey.fromUnit(r.unitString()), r.RepeatCount), r);
                });

        // now do downsampling
        downsampleSites();

        // now write only the downsampled sites to file
        String outputFile = RefGenomeMicrosatelliteFile.generateFilename(mConfig.OutputDir, mConfig.RefGenVersion);
        try (RefGenomeMicrosatelliteFile refGenomeMicrosatelliteFile = new RefGenomeMicrosatelliteFile(outputFile))
        {
            List<MicrosatelliteSite> sortedMicrosatelliteSites = new ArrayList<>(mDownSampledMicrosatelliteSites.values());
            sortedMicrosatelliteSites.sort(Comparator.comparing(site -> site.Region));
            sortedMicrosatelliteSites.forEach(refGenomeMicrosatelliteFile::writeRow);
        }

        RD_LOGGER.info("wrote {} microsatellite sites into {}", mDownSampledMicrosatelliteSites.size(), outputFile);

        Duration duration = Duration.between(start, Instant.now());

        RD_LOGGER.info("ref genome microsatellites finder took {}m {}s", duration.toMinutes(), duration.toSecondsPart());

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
            final Consumer<MicrosatelliteSite> refGenomeMsConsumer)
    {
        int chunkSize = 100_000;
        findMicrosatellites(referenceSequenceFile, minRepeatLength, refGenomeMsConsumer, chunkSize);
    }

    // algorithm to find short tandem repeats
    // at each base, we try find longest candidate starting from this base.
    // if a microsatellite is found, we start again from the base after.

    static void findMicrosatellites(
            final ReferenceSequenceFile referenceSequenceFile, int minNumRepeats,
            final Consumer<MicrosatelliteSite> refGenomeMsConsumer, int chunkSize)
    {
        MutableInt microsatelliteCounter = new MutableInt(0);

        List<SAMSequenceRecord> seqRecords = referenceSequenceFile.getSequenceDictionary().getSequences()
                .stream()
                .filter(o -> HumanChromosome.contains(o.getContig()))
                .collect(Collectors.toList());

        for(SAMSequenceRecord sequenceRecord : seqRecords)
        {
            RD_LOGGER.info("start processing chromosome {}", sequenceRecord.getContig());

            int length = sequenceRecord.getSequenceLength();

            // pending candidate, we do not accept a candidate when it is completed. We want to avoid
            // candidates that are too close to each other. They are dropped if too close.
            List<MicrosatelliteSite> pendingMicrosatellies = new ArrayList<>();

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
                            MicrosatelliteSite microsatelliteSite = new MicrosatelliteSite(sequenceRecord.getContig(),
                                    bestCandidate.startIndex,
                                    bestCandidate.startIndex + baseLength - 1, // change to inclusive
                                    bestCandidate.pattern);
                            pendingMicrosatellies.add(microsatelliteSite);

                            // check the panding microsatellites, see if any can be accepted
                            checkPendingMicrosatellites(pendingMicrosatellies, refGenomeMsConsumer, microsatelliteCounter);
                        }

                        bestCandidate = null;
                    }

                    // also start new Candidates at this location
                    // all these candidates have the current base as the last base of the pattern
                    for(int j = Math.max(i - JitterConstants.MAX_MICROSAT_UNIT_LENGTH + 1, 0); j <= i; ++j)
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

            RD_LOGGER.info("finished chromosome {}", sequenceRecord.getSequenceName());
        }

        RD_LOGGER.info("found {} microsatellite regions in ref genome", microsatelliteCounter);
    }

    // the aim of this code is to remove microsatellites that are too close to each other
    // We do not try to merge them for now, even though tools such as MsDetector would.
    // i.e. AAAAATAAAAAAA
    static void checkPendingMicrosatellites(
            final List<MicrosatelliteSite> pendingMicrosatellies, final Consumer<MicrosatelliteSite> refGenomeMsConsumer,
            final MutableInt microsatelliteCounter)
    {
        int groupStart = 0;

        // check the panding microsatellites, see if any can be accepted
        for(int i = 0; i < pendingMicrosatellies.size() - 1; ++i)
        {
            MicrosatelliteSite ms1 = pendingMicrosatellies.get(i);
            MicrosatelliteSite ms2 = pendingMicrosatellies.get(i + 1);

            if((ms2.referenceStart() - ms1.referenceEnd()) > JitterConstants.MIN_ADJACENT_MICROSAT_DISTANCE)
            {
                // previous group finished, if previous group only has 1 ms, we accept it, otherwise
                // remove them all from the pending list
                if(i == groupStart)
                {
                    // only 1 item, accept this
                    refGenomeMsConsumer.accept(ms1);
                    microsatelliteCounter.increment();
                    RD_LOGGER.trace("microsatellite: {}", ms1);
                }
                else
                {
                    // ms are too close to each other, remove them
                    for(int j = groupStart; j <= i; ++j)
                    {
                        RD_LOGGER.trace("reject microsatellite as too close to neighbour: {}", pendingMicrosatellies.get(j));
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
            final MicrosatelliteSite microsatelliteSite, final Map<String,List<GCProfile>> gcProfileMap)
    {
        List<GCProfile> gcProfiles = gcProfileMap.get(microsatelliteSite.chromosome());
        int siteMid = microsatelliteSite.Region.start() + (microsatelliteSite.Region.baseLength()) / 2;

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
            Validate.isTrue(gcProfile.overlaps(microsatelliteSite.Region.genomeRegion()));
            microsatelliteSite.setMappability(gcProfile.mappablePercentage());
        }
        else
        {
            microsatelliteSite.setMappability(0);
            RD_LOGGER.warn("microsatellite site({}) gc profile not found", microsatelliteSite.Region);
        }
    }

    public static int downsampleCount(UnitRepeatKey unitRepeatKey)
    {
        int rawValue = 25_000_000 / (int) Math.pow(2, (unitRepeatKey.NumRepeats + unitRepeatKey.Key.length * 2));
        return Math.max(MIN_TARGET_SITE_COUNT, rawValue);
    }

    // filter the microsatellites such that each type of (unit, length) is approximately the target count
    void downsampleSites() throws ExecutionException, InterruptedException
    {

        Map<String,List<BaseRegion>> regionsToKeep = ChrBaseRegion.loadChrBaseRegions(mConfig.BedFile);


        int numDigits = Integer.toString(mConfig.Threads - 1).length();
        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("thread-%0" + numDigits + "d").build();
        ExecutorService executorService = Executors.newFixedThreadPool(mConfig.Threads, namedThreadFactory);

        final List<CompletableFuture<Void>> futures = new ArrayList<>();

        // first we need to decide how much to downsample by. We do so by working out how many sites need to be kept
        for(UnitRepeatKey unitRepeatKey : mAllMicrosatelliteSites.keySet())
        {
            int downsampleCountPerType = downsampleCount(unitRepeatKey);
            RD_LOGGER.info("[{}] filtering microsatellite sites, target count per type = {}", unitRepeatKey, downsampleCountPerType);
            Runnable task = () -> downSampleMsType(downsampleCountPerType, regionsToKeep, unitRepeatKey);
            futures.add(CompletableFuture.runAsync(task, executorService));
        }

        // wait for completion
        CompletableFuture.allOf(futures.toArray(CompletableFuture[]::new)).get();

        RD_LOGGER.info("filtered {} microsatellite sites down to {}",
                mAllMicrosatelliteSites.size(), mDownSampledMicrosatelliteSites.size());
    }

    private static List<MicrosatelliteSite> downsampleList(List<MicrosatelliteSite> originalList, int numItemsToKeep, Random seed)
    {
        if(numItemsToKeep >= originalList.size())
            return originalList;
        else if(numItemsToKeep <= 0)
            return Collections.emptyList();

        Collections.shuffle(originalList, seed);
        return originalList.subList(0, numItemsToKeep);
    }

    private void downSampleMsType(
            final int targetCountPerType, final Map<String,List<BaseRegion>> regionsToKeep,
            final UnitRepeatKey unitRepeatKey)
    {
        Collection<MicrosatelliteSite> allList = mAllMicrosatelliteSites.get(unitRepeatKey);

        RD_LOGGER.info("[{}] filtering {} microsatellite sites", unitRepeatKey, allList.size());

        List<MicrosatelliteSite> sitesOutsideBedRegions = new ArrayList<>();
        List<MicrosatelliteSite> sitesInsideBedRegions = new ArrayList<>();

        // we filter to sites satisfying min mappability constraints
        // we then downsample to our target downsample count, prioritising sites inside the bed regions + 950 base buffer
        // if sites inside our bed regions already exceed target downsample count, we downsample these and take nothing from off-target

        int bedRegionExpansion = BED_REGION_EXPANSION;
        double mappabilityCutoff = MIN_MAPPABILITY;

        for(MicrosatelliteSite r : allList)
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
                sitesInsideBedRegions.add(r);
            }
            else
            {
                sitesOutsideBedRegions.add(r);
            }
        }

        // now decide how many to filter out from the rest
        int numSitesInBed = sitesInsideBedRegions.size();
        Random randomSeed = new Random(0);
        Collection<MicrosatelliteSite> filteredList = mDownSampledMicrosatelliteSites.get(unitRepeatKey);

        int sitesOutsideBedToKeep = targetCountPerType - numSitesInBed;
        List<MicrosatelliteSite> downsampledSitesInsideBedRegions = downsampleList(sitesInsideBedRegions, targetCountPerType, randomSeed);
        List<MicrosatelliteSite> downsampledSitesOutsideBedRegions = downsampleList(sitesOutsideBedRegions, sitesOutsideBedToKeep, randomSeed);
        filteredList.addAll(downsampledSitesInsideBedRegions);
        filteredList.addAll(downsampledSitesOutsideBedRegions);

        RD_LOGGER.info("[{}] filtered {} microsatellite sites down to {}, sites in bed: {} filtered to {}, sites outside bed: {} filtered to {}",
                unitRepeatKey, allList.size(), filteredList.size(), numSitesInBed, downsampledSitesInsideBedRegions.size(), sitesOutsideBedRegions.size(), downsampledSitesOutsideBedRegions.size());
    }

    public static void main(final String... args) throws Exception
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        Config.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        System.exit(new MicrosatelliteSiteFinder(configBuilder).run());
    }
}
