package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.primaryRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.secondMateRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.supplementaryConfig;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.supplementaryRecord;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.common.ContigEntry;
import com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TestGenome;
import com.hartwig.hmftools.tars.liftback.features.GenomicAlignmentScorer;
import com.hartwig.hmftools.tars.liftback.features.OverhangGate;
import com.hartwig.hmftools.tars.liftback.features.SupplementaryResolver;

import htsjdk.samtools.SAMRecord;

// Whole-pipeline test harness: declare the contig geometry, reference genome, annotated junctions and reads, then run() drives the
// real LiftBackGroupProcessor.processNameGroup (candidate features -> discriminator -> supplementary resolve -> mate patch -> NH /
// unmap policy) and Result asserts the lifted placement of every emitted record.
//
// Usage:
//   LiftBackScenario.create()
//       .contig(TarsTestFixtures.threeExonContig())
//       .genome(new TestGenome().with(CHR_1, 5000, 'A'))
//       .annotatedIntron(CHR_1, 200, 299)
//       .read(primary("frag1", TX_CONTIG, 51, "100M").mapQuality(0).as(151).xa("chr1_tx,+8725640,100M,0;"))
//       .read(mate("frag1", TX_CONTIG, 197, "50M"))
//       .read(supp("frag1", TX_CONTIG, 300, "50S50M").sa("chr1_tx,51,+,100M,60,0;"))
//       .run()
//       .assertLifted("frag1", ReadRole.PRIMARY, CHR_1, 150, "50M100N50M")
//       .assertSuppCount("frag1", 0);
public final class LiftBackScenario
{
    private final List<ContigEntry> mContigs = new ArrayList<>();
    private final Set<ChrBaseRegion> mAnnotatedIntrons = new HashSet<>();
    private final List<ReadSpec> mReads = new ArrayList<>();
    private TestGenome mGenome;
    private EnsemblAnnotationIndex mEnsemblAnnotationIndex;

    public enum ReadRole
    {
        PRIMARY,        // first-of-pair primary
        MATE,           // second-of-pair primary
        SUPPLEMENTARY   // a supplementary (first-of-pair)
    }

    public static LiftBackScenario create()
    {
        return new LiftBackScenario();
    }

    public LiftBackScenario contig(final ContigEntry entry)
    {
        mContigs.add(entry);
        return this;
    }

    public LiftBackScenario genome(final TestGenome genome)
    {
        mGenome = genome;
        return this;
    }

    public LiftBackScenario annotatedIntron(final String chromosome, final int start, final int end)
    {
        mAnnotatedIntrons.add(new ChrBaseRegion(chromosome, start, end));
        return this;
    }

    // Optional: drives the hidden-tie exon override in the MAPQ policy.
    public LiftBackScenario exonIndex(final EnsemblAnnotationIndex annotationIndex)
    {
        mEnsemblAnnotationIndex = annotationIndex;
        return this;
    }

    public LiftBackScenario read(final ReadSpec spec)
    {
        mReads.add(spec);
        return this;
    }

    public static ReadSpec primary(final String readName, final String contig, final int pos, final String cigar)
    {
        return new ReadSpec(readName, ReadRole.PRIMARY, primaryRecord(readName, contig, pos, cigar));
    }

    public static ReadSpec mate(final String readName, final String contig, final int pos, final String cigar)
    {
        return new ReadSpec(readName, ReadRole.MATE, secondMateRecord(readName, contig, pos, cigar));
    }

    public static ReadSpec supp(final String readName, final String contig, final int pos, final String cigar)
    {
        // the SA tag is left null; set it with .sa(...) when needed.
        return new ReadSpec(readName, ReadRole.SUPPLEMENTARY, supplementaryRecord(readName, contig, pos, cigar));
    }

    public Result run()
    {
        RefGenomeInterface ref = mGenome != null ? mGenome.asRefGenome() : null;
        EnsemblAnnotationIndex annotationIndex = EnsemblAnnotationIndex.fromJunctions(mAnnotatedIntrons);

        LiftBackDiscriminator resolver = mEnsemblAnnotationIndex != null
                ? new LiftBackDiscriminator(mContigs, mEnsemblAnnotationIndex)
                : new LiftBackDiscriminator(mContigs);
        SupplementaryResolver supplementary = new SupplementaryResolver(annotationIndex, ref, supplementaryConfig());
        OverhangGate overhangGate = new OverhangGate(ref);
        GenomicAlignmentScorer alignmentScorer = new GenomicAlignmentScorer(ref);

        LiftBackGroupProcessor processor = new LiftBackGroupProcessor(
                resolver, supplementary, overhangGate, alignmentScorer, ref, null);

        List<SAMRecord> emitted = new ArrayList<>();
        for(List<SAMRecord> group : groupByReadName())
        {
            processor.processNameGroup(group, emitted::add);
        }
        return new Result(emitted);
    }

    // reads sharing a name form one name-group, in declaration order.
    private List<List<SAMRecord>> groupByReadName()
    {
        List<List<SAMRecord>> groups = new ArrayList<>();
        List<String> order = new ArrayList<>();
        for(ReadSpec spec : mReads)
        {
            int index = order.indexOf(spec.ReadName);
            if(index < 0)
            {
                order.add(spec.ReadName);
                groups.add(new ArrayList<>());
                index = groups.size() - 1;
            }
            groups.get(index).add(spec.Record);
        }
        return groups;
    }

    // set bases(...) when a ref-dependent pass (overhang gate, supplementary motif scan) must compare read against genome.
    public static final class ReadSpec
    {
        final String ReadName;
        final ReadRole Role;
        final SAMRecord Record;

        ReadSpec(final String readName, final ReadRole role, final SAMRecord record)
        {
            ReadName = readName;
            Role = role;
            Record = record;
        }

        public ReadSpec mapQuality(final int mapQuality)
        {
            Record.setMappingQuality(mapQuality);
            return this;
        }

        public ReadSpec bases(final String sequence)
        {
            Record.setReadBases(TarsTestFixtures.bases(sequence));
            return this;
        }

        public ReadSpec xa(final String xa)
        {
            Record.setAttribute("XA", xa);
            return this;
        }

        public ReadSpec sa(final String sa)
        {
            Record.setAttribute("SA", sa);
            return this;
        }

        public ReadSpec as(final int alignmentScore)
        {
            Record.setAttribute("AS", alignmentScore);
            return this;
        }

        public ReadSpec nm(final int numMismatches)
        {
            Record.setAttribute("NM", numMismatches);
            return this;
        }

        public ReadSpec tag(final String name, final Object value)
        {
            Record.setAttribute(name, value);
            return this;
        }
    }

    // Emitted records keyed by (readName, role), with assertion helpers over the lifted placement.
    public static final class Result
    {
        private final List<SAMRecord> mEmitted;

        Result(final List<SAMRecord> emitted)
        {
            mEmitted = emitted;
        }

        public List<SAMRecord> emitted()
        {
            return mEmitted;
        }

        public SAMRecord record(final String readName, final ReadRole role)
        {
            for(SAMRecord record : mEmitted)
            {
                if(!record.getReadName().equals(readName))
                {
                    continue;
                }
                if(matchesRole(record, role))
                {
                    return record;
                }
            }
            return null;
        }

        private static boolean matchesRole(final SAMRecord record, final ReadRole role)
        {
            return switch(role)
            {
                case PRIMARY -> !record.getSupplementaryAlignmentFlag()
                        && (!record.getReadPairedFlag() || record.getFirstOfPairFlag());
                case MATE -> !record.getSupplementaryAlignmentFlag()
                        && record.getReadPairedFlag() && !record.getFirstOfPairFlag();
                case SUPPLEMENTARY -> record.getSupplementaryAlignmentFlag();
            };
        }

        public Result assertLifted(final String readName, final ReadRole role, final String chrom, final int pos, final String cigar)
        {
            SAMRecord record = record(readName, role);
            assertNotNull("no emitted record for " + readName + "/" + role, record);
            assertEquals(readName + "/" + role + " chrom", chrom, record.getReferenceName());
            assertEquals(readName + "/" + role + " pos", pos, record.getAlignmentStart());
            assertEquals(readName + "/" + role + " cigar", cigar, record.getCigarString());
            return this;
        }

        public Result assertMapQuality(final String readName, final ReadRole role, final int mapQuality)
        {
            SAMRecord record = record(readName, role);
            assertNotNull("no emitted record for " + readName + "/" + role, record);
            assertEquals(readName + "/" + role + " mapQuality", mapQuality, record.getMappingQuality());
            return this;
        }

        public Result assertNoXa(final String readName, final ReadRole role)
        {
            SAMRecord record = record(readName, role);
            assertNotNull("no emitted record for " + readName + "/" + role, record);
            assertNull(readName + "/" + role + " XA should be absent", record.getStringAttribute("XA"));
            return this;
        }

        public Result assertXa(final String readName, final ReadRole role, final String xa)
        {
            SAMRecord record = record(readName, role);
            assertNotNull("no emitted record for " + readName + "/" + role, record);
            assertEquals(readName + "/" + role + " XA", xa, record.getStringAttribute("XA"));
            return this;
        }

        public Result assertUnmapped(final String readName, final ReadRole role)
        {
            SAMRecord record = record(readName, role);
            assertNotNull("no emitted record for " + readName + "/" + role, record);
            assertEquals(readName + "/" + role + " should be unmapped", true, record.getReadUnmappedFlag());
            return this;
        }

        // mate cross-pointer, set by LiftBackRecordOps.patchMateFields against the mate's lifted primary.
        public Result assertMate(final String readName, final ReadRole role, final String mateChrom, final int matePos)
        {
            SAMRecord record = record(readName, role);
            assertNotNull("no emitted record for " + readName + "/" + role, record);
            assertEquals(readName + "/" + role + " mate chrom", mateChrom, record.getMateReferenceName());
            assertEquals(readName + "/" + role + " mate pos", matePos, record.getMateAlignmentStart());
            return this;
        }

        public Result assertMateUnmapped(final String readName, final ReadRole role)
        {
            SAMRecord record = record(readName, role);
            assertNotNull("no emitted record for " + readName + "/" + role, record);
            assertEquals(readName + "/" + role + " mate should be unmapped", true, record.getMateUnmappedFlag());
            return this;
        }

        public Result assertProperPair(final String readName, final ReadRole role, final boolean expected)
        {
            SAMRecord record = record(readName, role);
            assertNotNull("no emitted record for " + readName + "/" + role, record);
            assertEquals(readName + "/" + role + " proper-pair", expected, record.getProperPairFlag());
            return this;
        }

        // SA must be present, start at the given chromosome, and carry no _tx contig name.
        public Result assertSaGenomic(final String readName, final ReadRole role, final String chrom)
        {
            SAMRecord record = record(readName, role);
            assertNotNull("no emitted record for " + readName + "/" + role, record);
            String sa = record.getStringAttribute("SA");
            assertNotNull(readName + "/" + role + " SA should be present", sa);
            assertTrue(readName + "/" + role + " SA should start at " + chrom + ": " + sa, sa.startsWith(chrom + ","));
            assertTrue(readName + "/" + role + " SA should not reference a _tx contig: " + sa, !sa.contains("_tx"));
            return this;
        }

        public Result assertSuppCount(final String readName, final int expected)
        {
            int count = 0;
            for(SAMRecord record : mEmitted)
            {
                if(record.getReadName().equals(readName) && record.getSupplementaryAlignmentFlag())
                {
                    ++count;
                }
            }
            assertEquals(readName + " supplementary count", expected, count);
            return this;
        }

        public Result assertEmittedCount(final int expected)
        {
            assertEquals("total emitted records", expected, mEmitted.size());
            return this;
        }
    }
}
