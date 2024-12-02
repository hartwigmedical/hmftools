package com.hartwig.hmftools.bamtools.remapper;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.assembly.alignment.AlignData;
import com.hartwig.hmftools.esvee.assembly.alignment.Aligner;
import com.hartwig.hmftools.esvee.common.FileCommon;
import htsjdk.samtools.*;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.*;
import java.util.stream.Collectors;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static java.lang.String.format;

public class AltContigRemapper {
    private final AltContigRemapperConfig mConfig;

    public AltContigRemapper(final AltContigRemapperConfig config) {
        mConfig = config;
    }

    static boolean isRelevantDictionaryItem(SAMSequenceRecord samSequenceRecord) {
        return !samSequenceRecord.getSequenceName().toLowerCase().startsWith("hla");
    }

    static boolean hasAltReference(SAMRecord record) {
        return isHlaAltReference(record.getReferenceName());
    }

    static boolean mateHasAltReference(SAMRecord record) {
        return isHlaAltReference(record.getMateReferenceName());
    }

    static boolean isHlaAltReference(String referenceName) {
        return referenceName.toLowerCase().startsWith("hla-");
    }

    static boolean isSupplementaryAlignment(int flags) { return (flags & SAMFlag.SUPPLEMENTARY_ALIGNMENT.intValue()) != 0; }


    static SAMRecord createRemappedRecord(@NotNull final SAMRecord record, @NotNull final BwaMemAlignment alignment) {
        SAMRecord remappedRecord = record.deepCopy();
        remappedRecord.setReferenceIndex(alignment.getRefId());
        remappedRecord.setAlignmentStart(alignment.getRefStart());
        remappedRecord.setCigarString(alignment.getCigar());
        remappedRecord.setMappingQuality(alignment.getMapQual());
//        remappedRecord.setMateReferenceIndex(alignment.getMateRefId());
//        remappedRecord.setMateAlignmentStart(alignment.getMateRefStart());
        return remappedRecord;
    }

    public void run() {
        BT_LOGGER.info("starting alt contig remapper");
        long startTimeMs = System.currentTimeMillis();
        Aligner aligner = mConfig.aligner();
        Map<ImmutablePair<String, Integer>, ImmutablePair<Integer, Integer>> locationMap = new HashMap<>();
        List<SAMRecord> recordsNeedingMateRemapping = new ArrayList<>();

        try (SamReader samReader = SamReaderFactory.makeDefault().open(new File(mConfig.OrigBamFile))) {
            SAMFileHeader fileHeader = samReader.getFileHeader();
            BT_LOGGER.info("Header: " + fileHeader);
            SAMFileHeader newHeader = fileHeader.clone();

            SAMSequenceDictionary dictionaryWithHlaAltsRemoved = new SAMSequenceDictionary();
            fileHeader.getSequenceDictionary()
                    .getSequences()
                    .stream()
                    .filter(AltContigRemapper::isRelevantDictionaryItem)
                    .forEach(dictionaryWithHlaAltsRemoved::addSequence);
            newHeader.setSequenceDictionary(dictionaryWithHlaAltsRemoved);
            String interimOutputFileName = mConfig.OutputFile + ".unsorted";
            SAMFileWriter bamWriter = new SAMFileWriterFactory().makeBAMWriter(newHeader, false, new File(interimOutputFileName));

            samReader.forEach(record -> {
                if (hasAltReference(record)) {
                    if (record.isSecondaryOrSupplementary()) {
                        BT_LOGGER.info("Secondary alt contig remapper: " + record.getReadName());
                        // TODO
                    }
                    List<BwaMemAlignment> alignments = aligner.alignSequence(record.getReadBases());
                    if (alignments.isEmpty()) {
                          BT_LOGGER.info("no alignment found for " + record);
                          // TODO
                    } else {
                        // The original and new alignment locations need to be stored for use
                        // in adjusting records that have hla alt mate references.
                        // We only store the location for the first alignment returned.
                        BwaMemAlignment alignment0 = alignments.get(0);
                        ImmutablePair<String,Integer> originalLocation = new ImmutablePair<>(record.getReferenceName(), record.getAlignmentStart());
                        ImmutablePair<Integer,Integer> remappedLocation = new ImmutablePair<>(alignment0.getRefId(), alignment0.getRefStart());
                        locationMap.put(originalLocation, remappedLocation);

                        alignments.forEach(alignment -> {
                            SAMRecord remappedRecord = createRemappedRecord(record, alignment);
                            if (mateHasAltReference(remappedRecord)) {
                                recordsNeedingMateRemapping.add(remappedRecord);
                            } else {
                                bamWriter.addAlignment(remappedRecord);
                            }
                        });
                    }
                } else {
                    if (mateHasAltReference(record)) {
                        // Set it aside so that the mate reference can be adjusted
                        recordsNeedingMateRemapping.add(record);
                    } else {
                        // Write it out unaltered.
                        bamWriter.addAlignment(record);
                    }
                }
            });

            // Adjust those records that have hla alt references.
            recordsNeedingMateRemapping.forEach(record -> {
                ImmutablePair<String,Integer> originalLocation = new ImmutablePair<>(record.getMateReferenceName(), record.getMateAlignmentStart());
                ImmutablePair<Integer,Integer> remappedLocation = locationMap.get(originalLocation);
                if (remappedLocation != null) {
                    SAMRecord adjusted = record.deepCopy();
                    adjusted.setMateReferenceIndex(remappedLocation.getLeft());
                    adjusted.setMateAlignmentStart(remappedLocation.getRight());
                    bamWriter.addAlignment(adjusted);
                } else {
                    // TODO - error???
                }
            });

            // Write the records to file.
            bamWriter.close();

            // Sort the interim output.
            FileCommon.writeSortedBam(interimOutputFileName, mConfig.OutputFile, mConfig.BamToolPath, 1);

        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
        BT_LOGGER.info("Remapping complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    public static void main(@NotNull final String[] args) {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        AltContigRemapperConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        AltContigRemapperConfig config = new AltContigRemapperConfig(configBuilder);
        AltContigRemapper regionSlicer = new AltContigRemapper(config);
        regionSlicer.run();
    }
}
