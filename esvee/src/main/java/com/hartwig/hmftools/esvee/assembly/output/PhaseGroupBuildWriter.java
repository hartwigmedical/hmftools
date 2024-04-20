package com.hartwig.hmftools.esvee.assembly.output;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.StringJoiner;

import com.hartwig.hmftools.esvee.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;

public class PhaseGroupBuildWriter
{
    private final AssemblyConfig mConfig;

    private final BufferedWriter mWriter;

    public PhaseGroupBuildWriter(final AssemblyConfig config)
    {
        mConfig = config;
        mWriter = initialiseWriter();
    }

    public boolean enabled() { return mWriter != null; }
    public void close() { closeBufferedWriter(mWriter);}

    private BufferedWriter initialiseWriter()
    {
        if(!mConfig.WriteTypes.contains(WriteType.PHASE_GROUP_BUILDING))
            return null;

        try
        {
            BufferedWriter writer = createBufferedWriter(mConfig.outputFilename(WriteType.PHASE_GROUP_BUILDING));

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("Assembly");
            sj.add("LinkingAssembly");
            sj.add("LinkType");
            sj.add("PhaseGroupSize");

            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to initialise phase group writer: {}", e.toString());
            return null;
        }
    }

    public synchronized void writePhaseGroupBuild(
            final JunctionAssembly assembly, final JunctionAssembly linkingAssembly, final String linkType, int phaseGroupSize)
    {
        if(mWriter == null)
            return;

        try
        {
            String assemblyInfo = format("%s", assembly.junction().toString());
            String linkingAssemblyInfo = linkingAssembly != null ? format("%s", linkingAssembly.junction().toString()) : "";

            mWriter.write(format("%s\t%s\t%s\t%d", assemblyInfo, linkingAssemblyInfo, linkType, phaseGroupSize));
            mWriter.newLine();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write phase group building info: {}", e.toString());
        }
    }
}
