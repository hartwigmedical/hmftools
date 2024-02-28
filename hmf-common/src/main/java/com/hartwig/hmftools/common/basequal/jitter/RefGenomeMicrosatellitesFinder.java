package com.hartwig.hmftools.common.basequal.jitter;

import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.GC_PROFILE;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.GC_PROFILE_DESC;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.loadChrGcProfileMap;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import static htsjdk.samtools.util.SequenceUtil.N;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.gc.ImmutableGCProfile;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.commons.lang3.Validate;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.StringUtil;

// simple class to find microsatellites regions in reference genome
public class RefGenomeMicrosatellitesFinder
{
    public static final Logger sLogger = LogManager.getLogger(RefGenomeMicrosatellitesFinder.class);

    static class Candidate
    {
        int startIndex = 0;

        // where we are up to, it is also same as end
        int currentEndIndex = 0;
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

        byte nextBase()
        {
            return pattern[(currentEndIndex - startIndex) % pattern.length];
        }

        int length()
        {
            return currentEndIndex - startIndex;
        }

        int numFullUnits()
        {
            return length() / pattern.length;
        }
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

    public static void findMicrosatellites(ReferenceSequenceFile referenceSequenceFile, int minRepeatLength,
            Consumer<RefGenomeMicrosatellite> refGenomeMsConsumer)
    {
        int chunkSize = 100_000;
        findMicrosatellites(referenceSequenceFile, minRepeatLength, refGenomeMsConsumer, chunkSize);
    }

    // algorithm to find short tandem repeats
    // at each base, we try find longest candidate starting from this base.
    // if a microsatellite is found, we start again from the base after.
    //
    static void findMicrosatellites(ReferenceSequenceFile referenceSequenceFile, int minNumRepeats,
            Consumer<RefGenomeMicrosatellite> refGenomeMsConsumer, int chunkSize)
    {
        MutableInt microsatelliteCounter = new MutableInt(0);

        List<SAMSequenceRecord> seqRecords = referenceSequenceFile.getSequenceDictionary().getSequences()
                .stream()
                .filter(o -> HumanChromosome.contains(o.getContig()))
                .collect(Collectors.toList());

        for(SAMSequenceRecord sequenceRecord : seqRecords)
        {
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

            sLogger.info("finished chromosome {}", sequenceRecord.getSequenceName());
        }

        sLogger.info("found {} microsatellite regions in ref genome", microsatelliteCounter);
    }

    // the aim of this code is to remove microsatellites that are too close to each other
    // We do not try to merge them for now, even though tools such as MsDetector would.
    // i.e. AAAAATAAAAAAA
    static void checkPendingMicrosatellites(List<RefGenomeMicrosatellite> pendingMicrosatellies, Consumer<RefGenomeMicrosatellite> refGenomeMsConsumer,
            MutableInt microsatelliteCounter)
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
                    sLogger.trace("microsatellite: {}", ms1);
                }
                else
                {
                    // ms are too close to each other, remove them
                    for(int j = groupStart; j <= i; ++j)
                    {
                        sLogger.trace("reject microsatellite as too close to neighbour: {}", pendingMicrosatellies.get(j));
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

    // overload for easy testing
    static List<RefGenomeMicrosatellite> findMicrosatellites(ReferenceSequenceFile referenceSequenceFile, int minNumRepeats, int chunkSize)
    {
        List<RefGenomeMicrosatellite> refGenomeMicrosatellites = new ArrayList<>();
        findMicrosatellites(referenceSequenceFile, minNumRepeats, refGenomeMicrosatellites::add, chunkSize);
        return refGenomeMicrosatellites;
    }

    private static class Config
    {
        public final String refGenomeFile;
        public final String outputDir;

        public final RefGenomeVersion refGenomeVersion;

        public final String GcProfilePath;

        public Config(final ConfigBuilder configBuilder)
        {
            refGenomeFile = configBuilder.getValue(REF_GENOME);
            outputDir = parseOutputDir(configBuilder);
            refGenomeVersion = RefGenomeVersion.from(configBuilder);
            GcProfilePath = configBuilder.getValue(GC_PROFILE);
        }

        public static void registerConfig(final ConfigBuilder configBuilder)
        {
            addRefGenomeVersion(configBuilder);
            configBuilder.addPath(REF_GENOME, true, REF_GENOME_CFG_DESC + ", required when using CRAM files");
            configBuilder.addPath(GC_PROFILE, true, GC_PROFILE_DESC);

            addOutputDir(configBuilder);

            addLoggingOptions(configBuilder);

            addSpecificChromosomesRegionsConfig(configBuilder);
        }

        public boolean isValid()
        {
            checkCreateOutputDir(outputDir);
            return true;
        }
    }

    private static void findAndWriteRefGenomeMicrosatellites(Config config) throws Exception
    {
        Map<String,List<GCProfile>> gcProfiles = loadChrGcProfileMap(config.GcProfilePath);

        String outputFile = RefGenomeMicrosatelliteFile.generateFilename(config.outputDir, config.refGenomeVersion);

        // find all the polymers
        IndexedFastaSequenceFile refGenome = new IndexedFastaSequenceFile(new File(config.refGenomeFile));
        try (RefGenomeMicrosatelliteFile refGenomeMicrosatelliteFile = new RefGenomeMicrosatelliteFile(outputFile))
        {
            RefGenomeMicrosatellitesFinder.findMicrosatellites(refGenome, JitterAnalyserConstants.MIN_MICROSAT_UNIT_COUNT,
                    r -> {
                            populateMappability(r, gcProfiles);
                            refGenomeMicrosatelliteFile.writeRow(r);
                    });

            //filterSpecificRegions(refGenomeMicrosatellites);
        }

        sLogger.info("output written to {}", outputFile);
    }

    private static void populateMappability(RefGenomeMicrosatellite refGenomeMicrosatellite, Map<String,List<GCProfile>> gcProfileMap)
    {
        // populate the mappability
        List<GCProfile> gcProfiles = gcProfileMap.get(refGenomeMicrosatellite.chromosome());
        int siteMid = refGenomeMicrosatellite.genomeRegion.start() + (refGenomeMicrosatellite.genomeRegion.baseLength()) / 2;

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
            Validate.isTrue(gcProfile.overlaps(refGenomeMicrosatellite.genomeRegion.genomeRegion()));
            refGenomeMicrosatellite.mappability = gcProfile.mappablePercentage();
        }
        else
        {
            refGenomeMicrosatellite.mappability = 0.0;
            sLogger.warn("microsatellite site: {} gc profile not found.", refGenomeMicrosatellite.genomeRegion);
        }
    }

    public static void main(final String... args) throws Exception
    {
        ConfigBuilder configBuilder = new ConfigBuilder("ErrorProfile");
        Config.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        setLogLevel(configBuilder);

        Config config = new Config(configBuilder);

        findAndWriteRefGenomeMicrosatellites(config);

        System.exit(0);
    }
}
