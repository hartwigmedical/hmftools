package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarFromStr;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.NONE;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.TransExonRef;

import org.apache.commons.compress.utils.Lists;

public class ChimericReadCache
{
    private final IsofoxConfig mConfig;
    private final BufferedWriter mReadWriter;

    public ChimericReadCache(final IsofoxConfig config)
    {
        mConfig = config;
        mReadWriter = config.Fusions.WriteChimericReads || config.RunPerfChecks ? initialiseReadWriter() : null;
    }

    public void close()
    {
        closeBufferedWriter(mReadWriter);
    }

    private BufferedWriter initialiseReadWriter()
    {
        try
        {
            final String outputFileName = mConfig.formOutputFile("chimeric_reads.csv");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("ReadGroupCount,ReadId,FusionGroup,Chromosome,PosStart,PosEnd,Orientation,Cigar");
            writer.write(",Flags,HasSupplAlign,SuppData,BasesStart,BasesEnd,MateChr,MatePosStart");
            writer.write(",GeneSetStart,GeneSetEnd,GenicStart,GenicEnd,InterGeneSplit");
            writer.write(",MappedCoords,ScRegionsMatchedStart,ScRegionsMatchedEnd");
            writer.write(",TopTransMatch,TransExonData,UpperTopTransMatch,UpperTransExonData");
            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write chimeric read data: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeReadData(final String readId, final List<FusionRead> reads, final String groupStatus)
    {
        if(mReadWriter == null)
            return;

        try
        {
            for(final FusionRead read : reads)
            {
                mReadWriter.write(String.format("%d,%s,%s,%s,%d,%d,%d,%s",
                        reads.size(), readId, groupStatus, read.Chromosome,
                        read.posStart(), read.posEnd(), read.Orientation, read.Cigar));

                /*
                mReadWriter.write(String.format(",%s,%s,%s,%s,%s,%s,%d,%s,%d",
                        read.isFirstOfPair(), read.HasSuppAlignment, read.isReadReversed(), read.isProperPair(),
                        read.SuppData != null, read.ReadBases, read.flags(), read.MateChromosome, read.MatePosStart));
                */

                mReadWriter.write(String.format(",%d,%s,%s,%s,%s,%s,%d",
                        read.Flags, read.HasSuppAlignment, read.SuppData != null ? read.SuppData.asCsv() : "NONE",
                        read.BoundaryBases[SE_START], read.BoundaryBases[SE_END], read.MateChromosome, read.MatePosStart));

                mReadWriter.write(String.format(",%d,%d,%s,%s,%s",
                        read.GeneCollections[SE_START], read.GeneCollections[SE_END],
                        read.IsGenicRegion[SE_START], read.IsGenicRegion[SE_END], read.HasInterGeneSplit));

                StringJoiner coordsStr = new StringJoiner(ITEM_DELIM);

                for(int[] coord : read.MappedCoords)
                {
                    coordsStr.add(String.format("%d:%d", coord[SE_START], coord[SE_END]));
                }

                mReadWriter.write(String.format(",%s,%d,%d",
                        coordsStr.toString(), read.SoftClipLengths[SE_START], read.SoftClipLengths[SE_END]));

                // log the transcript exons affected, and the highest matching transcript
                StringJoiner transExonData = new StringJoiner(ITEM_DELIM);
                RegionMatchType topTransMatchType = read.getRegionMatchType(SE_START);

                if(topTransMatchType != NONE)
                {
                    for(final FusionTransExon transExonRef : read.getTransExonRefs(SE_START))
                    {
                        transExonData.add(String.format("%d:%d", transExonRef.TransId, transExonRef.ExonRank));
                    }
                }

                StringJoiner upperTransExonData = new StringJoiner(ITEM_DELIM);
                RegionMatchType upperTopTransMatchType = NONE;

                if(read.spansGeneCollections() && read.getTransExonRefs(SE_END) != null)
                {
                    upperTopTransMatchType = read.getRegionMatchType(SE_END);

                    for(FusionTransExon transExonRef : read.getTransExonRefs(SE_END))
                    {
                        upperTransExonData.add(String.format("%d:%d", transExonRef.TransId, transExonRef.ExonRank));
                    }
                }

                mReadWriter.write(String.format(",%s,%s,%s,%s",
                        topTransMatchType, transExonData.toString().isEmpty() ? "NONE" : transExonData.toString(),
                        upperTopTransMatchType, upperTransExonData.toString().isEmpty() ? "NONE" : upperTransExonData.toString()));

                mReadWriter.newLine();
            }

        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write chimeric read data: {}", e.toString());
            return;
        }
    }

    public static List<FusionReadGroup> loadChimericReads(final String inputFile)
    {
        final List<FusionReadGroup> readGroupList = Lists.newArrayList();

        // TODO - use FusionRead instead
        ISF_LOGGER.error("current unsupported");

        if(!Files.exists(Paths.get(inputFile)))
        {
            ISF_LOGGER.error("invalid chimeric reads file: {}", inputFile);
            return readGroupList;
        }

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(inputFile));

            // skip field names
            String line = fileReader.readLine();

            if (line == null)
            {
                ISF_LOGGER.error("empty chimeric reads file: {}", inputFile);
                return readGroupList;
            }

            final Map<String,Integer> fieldsMap  = createFieldsIndexMap(line, DELIMITER);

            int readId = fieldsMap.get("ReadId");
            int fusionGroup = fieldsMap.get("FusionGroup");
            int chr = fieldsMap.get("Chromosome");
            int posStart = fieldsMap.get("PosStart");
            int posEnd = fieldsMap.get("PosEnd");
            int cigar = fieldsMap.get("Cigar");
            int insertSize = fieldsMap.get("InsertSize");
            int flags = fieldsMap.get("Flags");
            int suppAlgn = fieldsMap.get("SuppAlign");
            int bases = fieldsMap.get("Bases");
            int mateChr = fieldsMap.get("MateChr");
            int matePosStart = fieldsMap.get("MatePosStart");
            int geneSetStart = fieldsMap.get("GeneSetStart");
            int geneSetEnd = fieldsMap.get("GeneSetEnd");
            int genicStart = fieldsMap.get("GenicStart");
            int genicEnd = fieldsMap.get("GenicEnd");
            int interGeneSplit = fieldsMap.get("InterGeneSplit");
            int mappedCoords = fieldsMap.get("MappedCoords");
            int scrmStart = fieldsMap.get("ScRegionsMatchedStart");
            int scrmEnd = fieldsMap.get("ScRegionsMatchedEnd");
            int topTransMatch = fieldsMap.get("TopTransMatch");
            int upperTopTransMatch = fieldsMap.get("UpperTopTransMatch");
            int transExonData = fieldsMap.get("TransExonData");
            int upperTransExonData = fieldsMap.get("UpperTransExonData");

            FusionReadGroup readGroup = null;

            while ((line = fileReader.readLine()) != null)
            {
                String[] items = line.split(DELIMITER, -1);

                try
                {
                    if(items[fusionGroup].contains("INVALID"))
                        continue;

                    ReadRecord read = new ReadRecord(
                            items[readId],
                            items[chr],
                            Integer.parseInt(items[posStart]),
                            Integer.parseInt(items[posEnd]),
                            items[bases],
                            cigarFromStr((items[cigar])),
                            Integer.parseInt(items[insertSize]),
                            Integer.parseInt(items[flags]),
                            items[mateChr],
                            Integer.parseInt(items[matePosStart]));

                    String saData = items[suppAlgn];

                    if(!saData.equals("NONE"))
                        read.setSuppAlignment(saData);

                    read.setGeneCollection(SE_START, Integer.parseInt(items[geneSetStart]), Boolean.parseBoolean(items[genicStart]));
                    read.setGeneCollection(SE_END, Integer.parseInt(items[geneSetEnd]), Boolean.parseBoolean(items[genicEnd]));

                    if(Boolean.parseBoolean(items[interGeneSplit]))
                        read.setHasInterGeneSplit();

                    read.getMappedRegionCoords().clear();
                    read.getMappedRegionCoords().addAll(parseMappedCoords(items[mappedCoords]));

                    read.getSoftClipRegionsMatched()[SE_START] = Integer.parseInt(items[scrmStart]);
                    read.getSoftClipRegionsMatched()[SE_END] = Integer.parseInt(items[scrmEnd]);

                    RegionMatchType matchType = RegionMatchType.valueOf(items[topTransMatch]);

                    if(matchType != NONE)
                    {
                        List<TransExonRef> transExonRefs = parseTransExonRefs(items[transExonData]);
                        // TODO - set into the map instead
                        // read.getReadTransExonRefs().put(matchType, transExonRefs);
                    }

                    matchType = RegionMatchType.valueOf(items[upperTopTransMatch]);

                    if(matchType != NONE)
                    {
                        List<TransExonRef> transExonRefs = parseTransExonRefs(items[upperTransExonData]);
                        // read.getReadTransExonRefs(SE_END).put(matchType, transExonRefs);
                    }

                    // reads are written by fragment so all reads will be sequential
                    if(readGroup == null || !readGroup.ReadId.equals(read.Id))
                    {
                        // readGroup = new FusionReadGroup(read);
                    }
                    else
                    {
                        // readGroup.Reads.add(read);

                        if(readGroup.isComplete())
                            readGroupList.add(readGroup);
                    }
                }
                catch (Exception e)
                {
                    ISF_LOGGER.error("failed to parse chimeric read data: {}", line);
                    return Lists.newArrayList();
                }
            }

            ISF_LOGGER.info("loaded {} chimeric fragment reads from file({})", readGroupList.size(), inputFile);
        }
        catch (IOException e)
        {
            ISF_LOGGER.warn("failed to load chimeric reads file({}): {}", inputFile, e.toString());
        }

        return readGroupList;
    }

    private static List<int[]> parseMappedCoords(final String data)
    {
        List<int[]> mappedCoords = Lists.newArrayList();

        for(String ref : data.split(ITEM_DELIM, -1))
        {
            String[] items = ref.split(":");
            if(items.length != 2)
                continue;

            mappedCoords.add(new int[] { Integer.parseInt(items[0]), Integer.parseInt(items[1]) });
        }

        return mappedCoords;
    }

    private static List<TransExonRef> parseTransExonRefs(final String data)
    {
        List<TransExonRef> transExonRefs = Lists.newArrayList();

        for(String ref : data.split(ITEM_DELIM, -1))
        {
            String[] items = ref.split(":");
            if(items.length != 4)
                continue;

            transExonRefs.add(new TransExonRef(items[0], Integer.parseInt(items[1]), items[2], Integer.parseInt(items[3])));
        }

        return transExonRefs;
    }

}
