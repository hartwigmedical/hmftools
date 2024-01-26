package com.hartwig.hmftools.esvee.output;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.WriteType;
import com.hartwig.hmftools.esvee.common.AssemblySupport;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.read.ReadUtils;

public class AssemblyWriter
{
    private final SvConfig mConfig;

    private final BufferedWriter mWriter;

    // write info about assemblies
    public AssemblyWriter(final SvConfig config)
    {
        mConfig = config;

        mWriter = initialiseWriter();
    }

    public void close() { closeBufferedWriter(mWriter);}

    private BufferedWriter initialiseWriter()
    {
        if(!mConfig.WriteTypes.contains(WriteType.ASSEMBLIES))
            return null;

        try
        {
            BufferedWriter writer = createBufferedWriter(mConfig.outputFilename(WriteType.ASSEMBLIES));

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("Chromosome");
            sj.add("JunctionPosition");
            sj.add("JunctionOrientation");
            sj.add("RangeStart");
            sj.add("RangeEnd");
            sj.add("Length");
            sj.add("SupportCount");
            sj.add("SoftClipMismatches");
            sj.add("RefBaseMismatches");
            sj.add("JunctionSequence");

            sj.add("AvgNmCount");
            sj.add("AvgIndelLength");
            sj.add("AvgRefSideSoftClip");
            sj.add("AvgBaseQual");
            sj.add("AvgMapQual");
            sj.add("InitialReadId");

            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to initialise assembly writer: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeAssembly(final JunctionAssembly assembly)
    {
        if(mWriter == null)
            return;

        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(assembly.initialJunction().Chromosome);
            sj.add(String.valueOf(assembly.initialJunction().Position));
            sj.add(String.valueOf(assembly.initialJunction().Orientation));
            sj.add(String.valueOf(assembly.minAlignedPosition()));
            sj.add(String.valueOf(assembly.maxAlignedPosition()));
            sj.add(String.valueOf(assembly.length()));

            int refBaseMismatches = 0;
            int softClipBaseMismatches = 0;

            for(AssemblySupport support : assembly.support())
            {
                refBaseMismatches += support.referenceMismatches();
                softClipBaseMismatches += support.junctionMismatches();
            }

            sj.add(String.valueOf(assembly.supportCount()));
            sj.add(String.valueOf(softClipBaseMismatches));
            sj.add(String.valueOf(refBaseMismatches));

            sj.add(assembly.formSequence(5));

            // ref sequence stats purely for analysis
            ReadStats readStats = buildReadStats(assembly.support(), assembly.initialJunction().isForward());
            sj.add(statString(readStats.NmCountTotal, assembly.supportCount()));
            sj.add(statString(readStats.IndelLengthTotal, assembly.supportCount()));
            sj.add(statString(readStats.RefSideSoftClipLengthTotal, assembly.supportCount()));
            sj.add(statString(readStats.BaseQualTotal, assembly.supportCount()));
            sj.add(statString(readStats.MapQualTotal, assembly.supportCount()));

            sj.add(assembly.initialRead().getName());

            mWriter.write(sj.toString());
            mWriter.newLine();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write assembly: {}", e.toString());
        }
    }

    private static String statString(int count, double readCount)
    {
        double avgValue = count/readCount;
        return format("%d", round(avgValue));
    }

    private class ReadStats
    {
        public int NmCountTotal;
        public int IndelLengthTotal;
        public int RefSideSoftClipLengthTotal;
        public int BaseQualTotal;
        public int MapQualTotal;

        public ReadStats()
        {
            NmCountTotal = 0;
            IndelLengthTotal = 0;
            RefSideSoftClipLengthTotal = 0;
            BaseQualTotal = 0;
            MapQualTotal = 0;
        }
    }

    private ReadStats buildReadStats(final List<AssemblySupport> supportReads, boolean isForwardJunction)
    {
        ReadStats readStats = new ReadStats();

        for(AssemblySupport support : supportReads)
        {
            Read read = support.read();
            readStats.NmCountTotal += read.numberOfEvents();
            readStats.MapQualTotal += read.mappingQuality();
            readStats.BaseQualTotal += ReadUtils.avgBaseQuality(read);
            readStats.IndelLengthTotal += read.cigarElements().stream().filter(x -> x.getOperator().isIndel()).mapToInt(x -> x.getLength()).sum();

            if(isForwardJunction)
                readStats.RefSideSoftClipLengthTotal += read.alignmentStart() - read.unclippedStart();
            else
                readStats.RefSideSoftClipLengthTotal += read.unclippedEnd() - read.alignmentEnd();
        }

        return readStats;
    }
}
