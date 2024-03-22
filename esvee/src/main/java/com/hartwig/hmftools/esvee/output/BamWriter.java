package com.hartwig.hmftools.esvee.output;

import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;

import java.io.File;
import java.io.FileOutputStream;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.common.AssemblyLink;
import com.hartwig.hmftools.esvee.common.AssemblySupport;
import com.hartwig.hmftools.esvee.common.IndelCoords;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.PhaseSet;
import com.hartwig.hmftools.esvee.common.SupportType;
import com.hartwig.hmftools.esvee.read.Read;

import htsjdk.samtools.BAMStreamWriter;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamWriter
{
    // private final BAMStreamWriter mWriter;
    private final SAMFileWriter mWriter;
    private SAMFileHeader mHeader;
    private boolean mValid;

    public BamWriter(final SvConfig config)
    {
        mHeader = null;
        mValid = true;
        mWriter = initialiseBamWriter(config);
        // mWriter = initialiseBam(config);
    }

    public boolean isValid() { return mWriter != null && mValid; }

    private BAMStreamWriter initialiseBam(final SvConfig config)
    {
        if(!config.WriteTypes.contains(WriteType.ASSEMBLY_BAM))
            return null;

        mHeader = new SAMFileHeader();
        mHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        RefGenomeSource refGenomeFile = (RefGenomeSource)config.RefGenome;

        refGenomeFile.refGenomeFile().getSequenceDictionary().getSequences().stream()
                .map(seq -> new SAMSequenceRecord(seq.getSequenceName(), seq.getSequenceLength()))
                .forEach(mHeader::addSequence);

        String outputFilename = config.outputFilename(WriteType.ASSEMBLY_BAM);

        try
        {
            BAMStreamWriter writer = new BAMStreamWriter(
                    new FileOutputStream(outputFilename),
                    new FileOutputStream(outputFilename + ".bai"),
                    null, 0, mHeader);

            // could take these from the original BAM(s)
            String programId = "ESVEE";
            String readGroupId = "NONE";
            mHeader.addProgramRecord(new SAMProgramRecord(programId));
            mHeader.addReadGroup(new SAMReadGroupRecord(readGroupId));

            writer.writeHeader(mHeader);

            return writer;
        }
        catch(final Exception e)
        {
            SV_LOGGER.error("failure while initialising output BAM", e.toString());
            mValid = false;
            return null;
        }
    }

    private SAMFileWriter initialiseBamWriter(final SvConfig config)
    {
        SAMFileHeader fileHeader = new SAMFileHeader();

        fileHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        RefGenomeSource refGenomeFile = (RefGenomeSource)config.RefGenome;

        refGenomeFile.refGenomeFile().getSequenceDictionary().getSequences().stream()
                .map(seq -> new SAMSequenceRecord(seq.getSequenceName(), seq.getSequenceLength()))
                .forEach(fileHeader::addSequence);

        String programId = "ESVEE";
        String readGroupId = "NONE";
        fileHeader.addProgramRecord(new SAMProgramRecord(programId));
        fileHeader.addReadGroup(new SAMReadGroupRecord(readGroupId));

        String outputBam = config.outputFilename(WriteType.ASSEMBLY_BAM);
        return new SAMFileWriterFactory().makeBAMWriter(fileHeader, false, new File(outputBam));
    }

    public void close()
    {
        if(mWriter != null)
        {
            mWriter.close();
            // mWriter.finish(true);
        }
    }

    private static final String ASSEMBLY_ID_TAG = "AS";
    private static final String READ_TYPE_TAG = "RT";
    private static final String ORIG_CIGAR_TAG = "OC";
    private static final String SPLIT_READ_COUNT = "SR";
    private static final String DISC_READ_COUNT = "DP";
    private static final String SAMPLE_TYPE_REF = "REF";
    private static final String SAMPLE_TYPE_TUMOR = "TUMOR";

    private static void setReadTypeTag(final SAMRecord read, final String sampleType, final String supportType)
    {
        read.setAttribute(READ_TYPE_TAG, sampleType + "_" + supportType);
    }

    public void writeAssembly(final JunctionAssembly assembly)
    {
        if(mWriter == null || !mValid)
            return;

        // write the assembly itself and then all its contributing reads
        AssemblyLink assemblyLink = null;

        if(assembly.phaseGroup() != null)
        {
            PhaseSet phaseSet = assembly.phaseGroup().findPhaseSet(assembly);

            if(phaseSet != null)
                assemblyLink = phaseSet.findSplitLink(assembly);
        }

        String assemblyReadId = format("ASSEMBLY_%d_%s", assembly.id(), assembly.junction().coords());

        int splitFragmentCount = 0;
        int discordantFragmentCount = 0;
        int totalMapQuality = 0;
        Set<String> uniqueFragmentIds = Sets.newHashSet();

        for(AssemblySupport support : assembly.support())
        {
            Read read = support.read();
            SAMRecord record = read.bamRecord();

            record.setAttribute(ASSEMBLY_ID_TAG, assemblyReadId);

            setReadTypeTag(record, read.isReference() ? SAMPLE_TYPE_REF : SAMPLE_TYPE_TUMOR, support.type().toString());

            if(read.originalCigarString() != read.cigarString())
                record.setAttribute(ORIG_CIGAR_TAG, read.originalCigarString());

            mWriter.addAlignment(record);

            if(!uniqueFragmentIds.contains(record.getReadName()))
            {
                uniqueFragmentIds.add(record.getReadName());

                if(support.type().isSplitSupport())
                {
                    ++splitFragmentCount;
                    totalMapQuality += read.mappingQuality();
                }
                else if(support.type() == SupportType.DISCORDANT)
                {
                    ++discordantFragmentCount;
                }
            }
        }

        int averageMapQual = (int)round(totalMapQuality / (double)splitFragmentCount);

        boolean isLinkedIndel = assembly.indel() && assemblyLink != null;

        if(!isLinkedIndel || assembly.isForwardJunction())
        {
            SAMRecord assemblyRead = createAssemblyRead(
                    assembly, assemblyReadId, assemblyLink, splitFragmentCount, discordantFragmentCount, averageMapQual);

            // mWriter.writeAlignment(assemblyRead);
            mWriter.addAlignment(assemblyRead);
        }
    }

    private static String formSplitAssemblyCigar(final JunctionAssembly assembly)
    {
        if(assembly.isForwardJunction())
        {
            return format("%dM%dS", assembly.refBaseLength(), assembly.extensionLength());
        }
        else
        {
            return format("%dS%dM", assembly.extensionLength(), assembly.refBaseLength());
        }
    }

    private SAMRecord createAssemblyRead(
            final JunctionAssembly assembly, final String assemblyReadId, final AssemblyLink assemblyLink,
            final int splitFragmentCount, final int discordantFragmentCount, final int avgMapQuality)
    {
        JunctionAssembly linkedAssembly = assemblyLink != null ? assemblyLink.otherAssembly(assembly) : null;

        boolean isLinkedIndel = assembly.indel() && linkedAssembly != null;

        SAMRecord record = new SAMRecord(mHeader);

        record.setReadName(assemblyReadId);
        record.setReadBases(assembly.bases());
        record.setBaseQualities(assembly.baseQuals());

        record.setMappingQuality(avgMapQuality);
        record.setReferenceName(assembly.junction().Chromosome);

        String cigarString;
        int alignmentStart;

        int mateAlignmentStart = NO_POSITION;
        String mateCigar = NO_CIGAR;

        if(isLinkedIndel)
        {
            IndelCoords indelCoords = assembly.initialRead().indelCoords();
            int upperRefLength = linkedAssembly.refBaseLength();

            alignmentStart = assembly.minAlignedPosition();

            cigarString = format("%dM%d%s%dM",
                    assembly.refBaseLength(),
                    indelCoords.isDelete() ? CigarOperator.D.toString() : CigarOperator.D.toString(),
                    indelCoords.Length, upperRefLength);
        }
        else
        {
            alignmentStart = assembly.isForwardJunction() ? assembly.minAlignedPosition() : assembly.junction().Position;
            cigarString = formSplitAssemblyCigar(assembly);

            if(linkedAssembly != null)
            {
                mateCigar = formSplitAssemblyCigar(linkedAssembly);
                mateAlignmentStart = linkedAssembly.isForwardJunction() ? linkedAssembly.minAlignedPosition() : linkedAssembly.junction().Position;
            }
        }

        record.setAlignmentStart(alignmentStart);
        record.setCigarString(cigarString);

        record.setFlags(assembly.initialRead().getFlags());

        if(linkedAssembly != null)
        {
            record.setMateReferenceName(linkedAssembly.junction().Chromosome);
            record.setMateAlignmentStart(mateAlignmentStart);
            // record.setMateReferenceIndex(templateRead.getMateReferenceIndex());
            record.setAttribute(MATE_CIGAR_ATTRIBUTE, mateCigar);
        }

        record.setDuplicateReadFlag(false);

        record.setInferredInsertSize(0);
        record.setAttribute(NUM_MUTATONS_ATTRIBUTE, 0);

        // special tags
        record.setAttribute(ASSEMBLY_ID_TAG, assemblyReadId);
        record.setAttribute(READ_TYPE_TAG, "ASSEMBLY");
        record.setAttribute(SPLIT_READ_COUNT, splitFragmentCount);
        record.setAttribute(DISC_READ_COUNT, discordantFragmentCount);

        return record;
    }
}
