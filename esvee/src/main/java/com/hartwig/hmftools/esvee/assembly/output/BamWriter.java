package com.hartwig.hmftools.esvee.assembly.output;

import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.esvee.assembly.output.WriteType.ASSEMBLY_BAM;
import static com.hartwig.hmftools.esvee.common.FileCommon.writeSortedBam;

import java.io.File;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.esvee.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.common.IndelCoords;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.PhaseSet;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;
import com.hartwig.hmftools.esvee.assembly.read.Read;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;

public class BamWriter
{
    private final AssemblyConfig mConfig;
    private final SAMFileWriter mWriter;
    private SAMFileHeader mHeader;

    private boolean mValid;
    private String mUnsortedBam;

    public BamWriter(final AssemblyConfig config)
    {
        mConfig = config;
        mValid = true;
        mUnsortedBam = null;

        if(config.WriteTypes.contains(WriteType.ASSEMBLY_BAM))
            mWriter = initialiseBamWriter(config);
        else
            mWriter = null;
    }

    public boolean isValid() { return mWriter != null && mValid; }

    private SAMFileWriter initialiseBamWriter(final AssemblyConfig config)
    {
        mHeader = new SAMFileHeader();

        mHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        RefGenomeSource refGenomeFile = (RefGenomeSource)config.RefGenome;

        refGenomeFile.refGenomeFile().getSequenceDictionary().getSequences().stream()
                .map(seq -> new SAMSequenceRecord(seq.getSequenceName(), seq.getSequenceLength()))
                .forEach(mHeader::addSequence);

        String programId = "ESVEE";
        String readGroupId = "NONE";
        mHeader.addProgramRecord(new SAMProgramRecord(programId));
        mHeader.addReadGroup(new SAMReadGroupRecord(readGroupId));

        String sortedBam = config.outputFilename(ASSEMBLY_BAM);
        mUnsortedBam = sortedBam.replaceAll(ASSEMBLY_BAM.fileId(), "unsorted." + ASSEMBLY_BAM.fileId());
        return new SAMFileWriterFactory().makeBAMWriter(mHeader, false, new File(mUnsortedBam));
    }

    public void close()
    {
        if(mWriter != null)
        {
            mWriter.close();

            String sortedBam = mConfig.outputFilename(ASSEMBLY_BAM);
            writeSortedBam(mUnsortedBam, sortedBam, mConfig.BamToolPath, mConfig.Threads);
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
        PhaseSet phaseSet = assembly.phaseSet();
        AssemblyLink assemblyLink = phaseSet != null ? phaseSet.findSplitLink(assembly) : null;

        String assemblyReadId = format("ASSEMBLY_%d_%s", assembly.id(), assembly.junction().coords());

        int splitFragmentCount = 0;
        int discordantFragmentCount = 0;
        int totalMapQuality = 0;
        Set<String> uniqueFragmentIds = Sets.newHashSet();

        for(SupportRead support : assembly.support())
        {
            Read read = support.cachedRead();
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
            IndelCoords indelCoords = assembly.indelCoords();
            int upperRefLength = linkedAssembly.refBaseLength();

            alignmentStart = indelMinAlignedReadPosition(assembly);

            cigarString = format("%dM%d%s%dM",
                    assembly.refBaseLength(),
                    indelCoords.isDelete() ? CigarOperator.D.toString() : CigarOperator.D.toString(),
                    indelCoords.Length, upperRefLength);
        }
        else
        {
            alignmentStart = assembly.minAlignedPosition();
            cigarString = formSplitAssemblyCigar(assembly);

            if(linkedAssembly != null)
            {
                mateCigar = formSplitAssemblyCigar(linkedAssembly);
                mateAlignmentStart = linkedAssembly.minAlignedPosition();
            }
        }

        record.setAlignmentStart(alignmentStart);
        record.setCigarString(cigarString);

        record.setFlags(assembly.support().get(0).flags());

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

    private static int indelMinAlignedReadPosition(final JunctionAssembly assembly)
    {
        return assembly.support().stream().mapToInt(x -> x.alignmentStart()).min().orElse(0);
    }
}
