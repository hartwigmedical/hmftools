package com.hartwig.hmftools.esvee.alignment;

import static com.hartwig.hmftools.common.region.SpecificRegions.parseStandardFormat;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.output.AlignmentWriter.FLD_ASSEMLY_INFO;
import static com.hartwig.hmftools.esvee.assembly.output.AlignmentWriter.FLD_CIGAR;
import static com.hartwig.hmftools.esvee.assembly.output.AlignmentWriter.FLD_FLAGS;
import static com.hartwig.hmftools.esvee.assembly.output.AlignmentWriter.FLD_XA_TAG;
import static com.hartwig.hmftools.esvee.assembly.output.AlignmentWriter.FLD_MAP_QUAL;
import static com.hartwig.hmftools.esvee.assembly.output.AlignmentWriter.FLD_MD_TAG;
import static com.hartwig.hmftools.esvee.assembly.output.AlignmentWriter.FLD_NMATCHES;
import static com.hartwig.hmftools.esvee.assembly.output.AlignmentWriter.FLD_REF_LOCATION;
import static com.hartwig.hmftools.esvee.assembly.output.AlignmentWriter.FLD_SCORE;
import static com.hartwig.hmftools.esvee.assembly.output.AlignmentWriter.FLD_SEQUENCE_COORDS;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class AlignmentCache
{
    private final Map<String,List<AlignData>> mAssemblyAlignmentData;

    private static final String ALIGNMENT_FILE = "alignment_file";

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(ALIGNMENT_FILE, false, "Pre-generated assembly alignment results");
    }

    public static String filename(final ConfigBuilder configBuilder) { return configBuilder.getValue(ALIGNMENT_FILE); }

    public AlignmentCache(final String filename)
    {
        mAssemblyAlignmentData = Maps.newHashMap();

        if(filename != null)
            loadAlignmentData(filename);

    }

    public boolean enabled() { return !mAssemblyAlignmentData.isEmpty(); }

    public List<AlignData> findAssemblyAlignments(final String assemblyInfo)
    {
        List<AlignData> alignments = mAssemblyAlignmentData.get(assemblyInfo);
        return alignments != null ? alignments : Collections.emptyList();
    }

    private void loadAlignmentData(final String filename)
    {
        if(filename == null || filename.isEmpty())
            return;

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));

            String header = lines.get(0);
            lines.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

            int assInfoIndex = fieldsIndexMap.get(FLD_ASSEMLY_INFO);
            int refLocationIndex = fieldsIndexMap.get(FLD_REF_LOCATION);
            int seqCoordsIndex = fieldsIndexMap.get(FLD_SEQUENCE_COORDS);
            int scoreIndex = fieldsIndexMap.get(FLD_SCORE);
            int flagIndex = fieldsIndexMap.get(FLD_FLAGS);
            int mqIndex = fieldsIndexMap.get(FLD_MAP_QUAL);
            int cigarIndex = fieldsIndexMap.get(FLD_CIGAR);
            int nmatchIndex = fieldsIndexMap.get(FLD_NMATCHES);
            int xaTagIndex = fieldsIndexMap.get(FLD_XA_TAG);
            int mdTagIndex = fieldsIndexMap.get(FLD_MD_TAG);

            List<AlignData> assemblyAlignments = null;
            String currentAssemblies = "";

            for(String line : lines)
            {
                String[] values = line.split(TSV_DELIM, -1);

                String refLocationStr = values[refLocationIndex];

                if(refLocationStr.startsWith("-1"))
                    continue;

                try
                {
                    ChrBaseRegion refLocation = parseStandardFormat(values[refLocationIndex]);
                    String assemblyInfo = values[assInfoIndex];
                    String[] sequenceCoords = values[seqCoordsIndex].split("-");
                    int seqStart = Integer.parseInt(sequenceCoords[0]);
                    int seqEnd = Integer.parseInt(sequenceCoords[1]);

                    int score = Integer.parseInt(values[scoreIndex]);
                    int mapQual = Integer.parseInt(values[mqIndex]);
                    int flags = Integer.parseInt(values[flagIndex]);

                    String cigar = values[cigarIndex];
                    int nMatches = Integer.parseInt(values[nmatchIndex]);
                    String xaTag = values[xaTagIndex];
                    String mdTag = values[mdTagIndex];

                    if(!currentAssemblies.equals(assemblyInfo))
                    {
                        currentAssemblies = assemblyInfo;
                        assemblyAlignments = Lists.newArrayList();
                        mAssemblyAlignmentData.put(assemblyInfo, assemblyAlignments);
                    }

                    assemblyAlignments.add(new AlignData(refLocation, seqStart, seqEnd, mapQual, score, flags, cigar, nMatches, xaTag, mdTag));
                }
                catch(Exception e)
                {
                    SV_LOGGER.error("invalid alignment data: {}", line);
                }
            }

            SV_LOGGER.info("loaded {} assembly alignments from file: {}",
                    mAssemblyAlignmentData.values().stream().mapToInt(x -> x.size()).sum(), filename);

        }
        catch(IOException exception)
        {
            SV_LOGGER.error("failed to read alignment data file({}): {}", filename, exception.toString());
        }
    }
}
