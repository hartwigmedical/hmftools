package com.hartwig.hmftools.bamtools.remapper;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.immune.ImmuneRegions;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.common.FileCommon;

import htsjdk.samtools.*;

import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;

public class AltContigRemapper
{
    private final AltContigRemapperConfig mConfig;

    public AltContigRemapper(final AltContigRemapperConfig config)
    {
        mConfig = config;
    }

    static boolean isRelevantDictionaryItem(SAMSequenceRecord samSequenceRecord)
    {
        return !samSequenceRecord.getSequenceName().toLowerCase().startsWith("hla");
    }

    public void run()
    {
        Set<String> readsWithDifferences = runToGetNamesOfReadsWithDifferences();
        BT_LOGGER.info("starting alt contig remapper");
        long startTimeMs = System.currentTimeMillis();

        Map<String, SAMRecord> differencesAlts1 = new HashMap<>();
        Map<String, SAMRecord> differencesAlts2 = new HashMap<>();
        Map<String, SAMRecord> differencesNoAlts1 = new HashMap<>();
        Map<String, SAMRecord> differencesNoAlts2 = new HashMap<>();

        try(SamReader samReader = SamReaderFactory.makeDefault()
                .open(new File("/Users/timlavers/work/scratch/comparison/no_alts.final.bam")))
        {
            samReader.forEach(record ->
            {
                if(!record.isSecondaryOrSupplementary())
                {
                    if(differencesNoAlts1.containsKey(record.getReadName()))
                    {
                        differencesNoAlts2.put(record.getReadName(), record);
                    }
                    else
                    {
                        differencesNoAlts1.put(record.getReadName(), record);
                    }
                }
            });
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }

        try(SamReader samReader = SamReaderFactory.makeDefault().open(new File(mConfig.OrigBamFile)))
        {

            // The header in the rewritten file needs to be the same
            // as the header in the original file but with the hla dictionary items removed.
            SAMFileHeader fileHeader = samReader.getFileHeader();
            SAMFileHeader newHeader = fileHeader.clone();
            SAMSequenceDictionary dictionaryWithHlaAltsRemoved = new SAMSequenceDictionary();
            fileHeader.getSequenceDictionary()
                    .getSequences()
                    .stream()
                    .filter(AltContigRemapper::isRelevantDictionaryItem)
                    .forEach(dictionaryWithHlaAltsRemoved::addSequence);
            newHeader.setSequenceDictionary(dictionaryWithHlaAltsRemoved);

            final BwaHlaRecordAligner aligner = new BwaHlaRecordAligner(mConfig.aligner(), newHeader);
            HlaTransformer transformer = new HlaTransformer(aligner);
            samReader.forEach(record ->
            {
                if(!record.isSecondaryOrSupplementary())
                {
                    transformer.process(record).forEach(alignment ->
                    {
                        if(!alignment.isSecondaryOrSupplementary())
                        {
                            if(readsWithDifferences.contains(alignment.getReadName()))
                            {
                                if(differencesAlts2.containsKey(alignment.getReadName()))
                                {
                                    differencesAlts1.put(alignment.getReadName(), alignment);
                                }
                                else
                                {
                                    differencesAlts2.put(alignment.getReadName(), alignment);
                                }
                            }
                        }
                    });
                }
            });
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
        readsWithDifferences.stream().sorted().forEach(name ->
        {
            System.out.println("------------------------------------------------------------------");
            System.out.println(name);
            System.out.println("no alts");
            printPair(differencesNoAlts1.get(name), differencesNoAlts2.get(name));
            System.out.println();

            System.out.println("remapped");
            printPair(differencesAlts1.get(name), differencesAlts2.get(name));
            System.out.println("------------------------------------------------------------------");
            System.out.println();
        });
        BT_LOGGER.info("Remapping complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    public Set<String> runToGetNamesOfReadsWithDifferences()
    {
        final AtomicInteger numberCompared = new AtomicInteger();
        final AtomicInteger numberMatching = new AtomicInteger();
        final AtomicInteger hlaRegionRemapping = new AtomicInteger();
        BT_LOGGER.info("starting alt contig remapper");
        long startTimeMs = System.currentTimeMillis();

        Map<String, SAMRecord> noAltsNegativeReads = new HashMap<>();
        Map<String, SAMRecord> noAltsPositiveReads = new HashMap<>();
        Set<String> namesOfReadsWithDifferences = new HashSet<>();

        try(SamReader samReader = SamReaderFactory.makeDefault()
                .open(new File("/Users/timlavers/work/scratch/comparison/no_alts.final.bam")))
        {
            samReader.forEach(record ->
            {
                if(!record.isSecondaryOrSupplementary())
                {
                    if(record.getReadNegativeStrandFlag())
                    {
                        noAltsNegativeReads.put(record.getReadName(), record);
                    }
                    else
                    {
                        noAltsPositiveReads.put(record.getReadName(), record);
                    }
                }
            });
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }

        try(SamReader samReader = SamReaderFactory.makeDefault().open(new File(mConfig.OrigBamFile)))
        {

            // The header in the rewritten file needs to be the same
            // as the header in the original file but with the hla dictionary items removed.
            SAMFileHeader fileHeader = samReader.getFileHeader();
            SAMFileHeader newHeader = fileHeader.clone();
            SAMSequenceDictionary dictionaryWithHlaAltsRemoved = new SAMSequenceDictionary();
            fileHeader.getSequenceDictionary()
                    .getSequences()
                    .stream()
                    .filter(AltContigRemapper::isRelevantDictionaryItem)
                    .forEach(dictionaryWithHlaAltsRemoved::addSequence);
            newHeader.setSequenceDictionary(dictionaryWithHlaAltsRemoved);

            // The initial output is unsorted.
            String interimOutputFileName = mConfig.OutputFile + ".unsorted";
            File interimOutputFile = new File(interimOutputFileName);
            SAMFileWriter bamWriter = new SAMFileWriterFactory().makeBAMWriter(newHeader, false, interimOutputFile);

            final BwaHlaRecordAligner aligner = new BwaHlaRecordAligner(mConfig.aligner(), newHeader);
            HlaTransformer transformer = new HlaTransformer(aligner);
            samReader.forEach(record ->
            {
                boolean isHla = HlaTransformer.hasSomeHlaReference(record);
                transformer.process(record).forEach(alignment ->
                {
                    bamWriter.addAlignment(alignment);
                    if (isHla)
                    {

                        if(!record.isSecondaryOrSupplementary() && !alignment.isSecondaryOrSupplementary())
                        {
                            SAMRecord noAltsRecord = alignment.getReadNegativeStrandFlag()
                                    ? noAltsNegativeReads.get(record.getReadName())
                                    : noAltsPositiveReads.get(record.getReadName());
                            if(noAltsRecord != null)
                            {
                                int qualDiff = Math.abs(alignment.getFlags() - noAltsRecord.getFlags());
                                numberCompared.getAndIncrement();
                                if(
                                        alignment.getAlignmentStart() != noAltsRecord.getAlignmentStart() ||
                                                //                                                    alignment.getFlags() != noAltsRecord.getFlags() ||
                                                qualDiff > 0 ||
                                                !Objects.equals(alignment.getCigarString(), noAltsRecord.getCigarString()) ||
                                                !Objects.equals(alignment.getReferenceIndex(), noAltsRecord.getReferenceIndex())
                                    //                                                    !Objects.equals(alignment.getMappingQuality(), noAltsRecord.getMappingQuality())
                                    //                                                    !Objects.equals(alignment.getSAMString(), noAltsRecord.getSAMString())
                                )
                                {
                                    String hlaOriginal = hlaRegion(noAltsRecord.getAlignmentStart());
                                    String hlaNew = hlaRegion(alignment.getAlignmentStart());
                                    if(hlaOriginal != null && hlaNew != null && !hlaOriginal.equals(hlaNew))
                                    {
                                        hlaRegionRemapping.incrementAndGet();
                                        //                                        } else {
                                    }
                                    namesOfReadsWithDifferences.add(record.getReadName());
                                    //                                        }
                                }
                                else
                                {
                                    numberMatching.getAndIncrement();
                                }
                            }
                            else
                            {
                                System.out.println("Not found or bad: " + record.getReadName());
                            }
                        }
                    }
                });
            });

            System.out.println("number compared: " + numberCompared.get());
            System.out.println("number matching: " + numberMatching.get());
            System.out.println("number of differences due to HLA region change: " + hlaRegionRemapping.get());
            // Deal with any unmatched reads.
            // Don't map these - log an error and write them out as they are
            List<SAMRecord> unmatched = transformer.unmatchedRecords();
            if(!unmatched.isEmpty())
            {
                BT_LOGGER.warn("Some HLA contig records were unmatched. " + unmatched);
                unmatched.forEach(bamWriter::addAlignment);
            }
            else
            {
                BT_LOGGER.info("No HLA contig records were unmatched.");
            }

            // Write the records to file.
            bamWriter.close();
            BT_LOGGER.info("BAM Writer closed.");

            // If the samtools path has been provided, sort the output. Else simply rename the unsorted file.
            if(mConfig.BamToolPath != null)
            {
                BT_LOGGER.info("Output file is to be sorted...");
                FileCommon.writeSortedBam(interimOutputFileName, mConfig.OutputFile, mConfig.BamToolPath, 1);
                BT_LOGGER.info("Sorting complete.");
            }
            else
            {
                File outputFile = new File(mConfig.OutputFile);
                boolean renamed = interimOutputFile.renameTo(outputFile);
                if(!renamed)
                {
                    BT_LOGGER.warn("Could not rename " + interimOutputFile + " to " + outputFile);
                }
            }
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
        BT_LOGGER.info("Remapping complete, mins({})", runTimeMinsStr(startTimeMs));
        return namesOfReadsWithDifferences;
    }

    private static String hlaRegion(int position)
    {
        List<ChrBaseRegion> regions = ImmuneRegions.getHlaRegions(RefGenomeVersion.V38);
        if(regions.get(0).containsPosition(position))
        {
            return "HLA-A";
        }
        if(regions.get(1).containsPosition(position))
        {
            return "HLA-B";
        }
        if(regions.get(2).containsPosition(position))
        {
            return "HLA-C";
        }

        return null;
    }

    private static void printPair(SAMRecord sam1, SAMRecord sam2)
    {
        if(sam1 != null && sam1.getReadNegativeStrandFlag())
        {
            System.out.println(summary(sam2));
            System.out.println(summary(sam1));
        }
        else
        {
            System.out.println(summary(sam1));
            System.out.println(summary(sam2));
        }

    }
    private static String summary(SAMRecord record)
    {
        if(record == null)
        {
            return null;
        }
        return record.getReferenceIndex() + "\t"
                + record.getAlignmentStart() + "\t"
                + hlaRegion(record.getAlignmentStart()) + "\t"
                + record.getFlags() + "\t"
                + record.getMappingQuality() + "\t"
                + record.getMateAlignmentStart() + "\t"
                + record.getCigarString() + "\t";
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        AltContigRemapperConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        AltContigRemapperConfig config = new AltContigRemapperConfig(configBuilder);
        new AltContigRemapper(config).run();
    }
}
