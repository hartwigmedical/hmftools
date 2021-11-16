package com.hartwig.hmftools.isofox.fusion.cohort;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.fusion.FusionJunctionType.KNOWN;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;
import com.hartwig.hmftools.isofox.fusion.FusionData;

import org.apache.commons.compress.utils.Lists;

public class UnknownSpliceAnalyser
{
    private final List<ChrBaseRegion> mKnownLineElements;
    private final BufferedWriter mWriter;

    public UnknownSpliceAnalyser(final FusionCohortConfig config, final BufferedWriter writer)
    {
        mWriter = writer;
        mKnownLineElements = Lists.newArrayList();

        loadLineElementsFile(config.LineElementsFile);
    }

    public void compareFusions(final String sampleId, final List<FusionData> fusions)
    {
        for(FusionData fusionData : fusions)
        {
            if(fusionData.JunctionTypes[FS_DOWN] != KNOWN)
                continue;

            if(!fusionData.GeneNames[FS_UP].isEmpty())
                continue;

            ChrBaseRegion lineSite = findLineSiteMatch(fusionData.Chromosomes[FS_UP], fusionData.JunctionPositions[FS_UP]);

            if(lineSite == null)
                continue;

            try
            {
                mWriter.write(String.format("%s,%d,%s,%s,%s,%d",
                        sampleId, fusionData.Id, fusionData.GeneIds[FS_DOWN], fusionData.GeneNames[FS_DOWN],
                        fusionData.Chromosomes[FS_DOWN], fusionData.JunctionPositions[FS_DOWN]));

                mWriter.write(String.format(",%s,%d,%d,%d,%d",
                        fusionData.Chromosomes[FS_UP], fusionData.JunctionPositions[FS_UP],
                        lineSite.start(), lineSite.end(), fusionData.TotalFrags));

                mWriter.newLine();
            }
            catch(IOException e)
            {
                ISF_LOGGER.error("failed to write splice-to-unknown data: {}", e.toString());
            }
        }
    }

    public static BufferedWriter initialiseWriter(final CohortConfig config)
    {
        try
        {
            final String outputFileName = config.formCohortFilename("splice_to_unknown.csv");

            final BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("SampleId,FusionId,GeneId,GeneName,FusionChromosome,FusionPosition");
            writer.write(",OtherChromosome,OtherJuncPosition,LineStart,LineEnd,TotalFragments");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to initialiser splice-to-unknown file: {}", e.toString());
            return null;
        }
    }

    private static final int LINE_MATCH_BUFFER = 1000;

    private ChrBaseRegion findLineSiteMatch(final String chromosome, final int position)
    {
        return mKnownLineElements.stream()
                .filter(x -> x.Chromosome.equals(chromosome))
                .filter(x -> positionWithin(position, x.start() - LINE_MATCH_BUFFER, x.end() + LINE_MATCH_BUFFER))
                .findFirst().orElse(null);
    }

    private void loadLineElementsFile(final String filename)
    {
        if(filename == null || filename.isEmpty())
            return;

        try
        {
            List<String> fileContents = Files.readAllLines(new File(filename).toPath());
            String header = fileContents.get(0);
            fileContents.remove(0);

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, ",");

            for(final String line : fileContents)
            {
                String[] values = line.split(",");

                final ChrBaseRegion lineRegion = new ChrBaseRegion(
                        values[fieldsIndexMap.get("Chromosome")],
                        Integer.parseInt(values[fieldsIndexMap.get("PosStart")]),
                        Integer.parseInt(values[fieldsIndexMap.get("PosEnd")]));

                mKnownLineElements.add(lineRegion);
            }

            ISF_LOGGER.info("loaded {} known line elements from file: {}", mKnownLineElements.size(), filename);
        }
        catch(IOException exception)
        {
            ISF_LOGGER.error("failed to read line element CSV file({})", filename);
        }
    }

}
