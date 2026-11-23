package com.hartwig.hmftools.pave.impact;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_TRANS_NAME;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.Nullable;

public class SpecificSpliceRegion
{
    public final String GeneName;
    public final String TranscriptName;
    public final int ExonRank;
    public final boolean UseStart;

    // in base terms relative to the exon boundary
    // eg -30 to -10 would be the 21 bases from 30 bases before the exon boundary to 10 bases before
    public final int BaseStart;
    public final int BaseEnd;

    public SpecificSpliceRegion(
            final String geneName, final String transcriptName, final int exonRank, final boolean useStart,
            final int baseStart, final int baseEnd)
    {
        GeneName = geneName;
        TranscriptName = transcriptName;
        ExonRank = exonRank;
        UseStart = useStart;
        BaseStart = baseStart;
        BaseEnd = baseEnd;
    }

    private static final String SPECIAL_SPLICE_REGIONS_FILE = "/splice/special_splice_regions.tsv";

    public static List<SpecificSpliceRegion> load(@Nullable final String filename)
    {
        List<SpecificSpliceRegion> spliceRegions = Lists.newArrayList();

        try
        {
            List<String> lines = Lists.newArrayList();

            if(filename != null)
            {
                lines.addAll(Files.readAllLines(new File(filename).toPath()));
            }
            else
            {
                lines.addAll(new BufferedReader(new InputStreamReader(
                        SpecificSpliceRegion.class.getResourceAsStream(SPECIAL_SPLICE_REGIONS_FILE)))
                        .lines().collect(Collectors.toList()));
            }

            String header = lines.get(0);
            lines.remove(0);

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
            int geneIndex = fieldsIndexMap.get(FLD_GENE_NAME);
            int transIndex = fieldsIndexMap.get(FLD_TRANS_NAME);
            int exonIndex = fieldsIndexMap.get("ExonRank");
            int isStartIndex = fieldsIndexMap.get("UseStart");
            int baseStartIndex = fieldsIndexMap.get("BaseStart");
            int baseEndIndex = fieldsIndexMap.get("BaseEnd");

            for(String line : lines)
            {
                String[] values = line.split(TSV_DELIM, -1);

                SpecificSpliceRegion region = new SpecificSpliceRegion(
                        values[geneIndex], values[transIndex], Integer.parseInt(values[exonIndex]), Boolean.parseBoolean(values[isStartIndex]),
                        Integer.parseInt(values[baseStartIndex]), Integer.parseInt(values[baseEndIndex]));
                spliceRegions.add(region);
            }
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to read fragile site file: {}", e.toString());
        }

        return spliceRegions;
    }
}
