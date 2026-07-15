package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.primaryRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.secondMateRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.supplementaryRecord;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.common.ContigEntry;
import com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TestGenome;
import com.hartwig.hmftools.tars.liftback.overhang.OverhangGate;
import com.hartwig.hmftools.tars.liftback.supplementary.AnnotatedJunctionIndex;
import com.hartwig.hmftools.tars.liftback.supplementary.SupplementaryResolver;
import com.hartwig.hmftools.tars.liftback.supplementary.RefSequenceSource;
import com.hartwig.hmftools.tars.liftback.supplementary.SupplementaryConfig;

import htsjdk.samtools.SAMRecord;

// Parametric, standalone whole-pipeline test harness. Declare the transcript-contig geometry, the reference genome,
// the annotated junctions and a set of reads (primary + mate + any supplementaries, each with arbitrary tags), then
// run() drives the real LiftBackGroupProcessor.processNameGroup exactly as production does (overhang gate ->
// discriminator -> supplementary resolve -> reclaim -> mate patch -> NH / unmap policy). Assert the lifted
// placement of every emitted record.
//
// Usage:
//   LiftBackScenario.create()
//       .contig(TarsTestFixtures.threeExonContig())
//       .genome(new TestGenome().with(CHR_1, 5000, 'A'))
//       .annotatedIntron(CHR_1, 200, 299)
//       .read(primary("frag1", TX_CONTIG, 51, "100M").mapq(0).as(151).xa("chr1_tx,+8725640,100M,0;"))
//       .read(mate("frag1", TX_CONTIG, 197, "50M"))
//       .read(supp("frag1", TX_CONTIG, 300, "50S50M").sa("chr1_tx,51,+,100M,60,0;"))
//       .run()
//       .assertLifted("frag1", ReadRole.PRIMARY, CHR_1, 150, "50M100N50M")
//       .assertMapq("frag1", ReadRole.PRIMARY, 60)
//       .assertNoXa("frag1", ReadRole.PRIMARY)
//       .assertSuppCount("frag1", 0);
public final class LiftBackScenario
{
    private final List<ContigEntry> mContigs = new ArrayList<>();
    private final Set<ChrBaseRegion> mAnnotatedIntrons = new HashSet<>();
    private final List<ReadSpec> mReads = new ArrayList<>();
    private TestGenome mGenome;
    private ExonRegionIndex mExonIndex;

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

    // optional: drives the hidden-tie exon override in the MAPQ policy.
    public LiftBackScenario exonIndex(final ExonRegionIndex exonIndex)
    {
        mExonIndex = exonIndex;
        return this;
    }

    public LiftBackScenario read(final ReadSpec spec)
    {
        mReads.add(spec);
        return this;
    }

    // ----- read factories: build a record for the given role, then chain tag setters -----

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
        // supplementaryRecord(readName, contig, pos, cigar) leaves the SA tag null; set it with .sa(...) when needed.
        return new ReadSpec(readName, ReadRole.SUPPLEMENTARY, supplementaryRecord(readName, contig, pos, cigar));
    }

    // ----- run the full pipeline and collect every emitted record -----

    public Result run()
    {
        RefSequenceSource ref = mGenome != null ? mGenome.asRefSource() : null;
        AnnotatedJunctionIndex junctionIndex = new AnnotatedJunctionIndex(mAnnotatedIntrons);

        LiftBackResolver resolver = mExonIndex != null
                ? new LiftBackResolver(mContigs, mExonIndex)
                : new LiftBackResolver(mContigs);
        SupplementaryResolver supplementary = new SupplementaryResolver(junctionIndex, ref, SupplementaryConfig.enabledDefaults());
        OverhangGate overhangGate = new OverhangGate(ref);

        LiftBackGroupProcessor processor = new LiftBackGroupProcessor(
                resolver, supplementary, overhangGate, ref, null, new LiftBackStats());

        List<SAMRecord> emitted = new ArrayList<>();
        for(List<SAMRecord> group : groupByReadName())
        {
            processor.processNameGroup(group, new LiftedMateInfoCache(), (record, result) -> emitted.add(record));
        }
        return new Result(emitted);
    }

    // reads sharing a name form one name-group; preserve declaration order across groups.
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

    // A record plus its role, with fluent tag setters. SEQ defaults to the fixture bases; override with bases(...)
    // when a ref-dependent pass (overhang gate / supplementary-resolve ref-verify / canonicalize) must compare read vs genome.
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

        public ReadSpec mapq(final int mapq)
        {
            Record.setMappingQuality(mapq);
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

        public Result assertMapq(final String readName, final ReadRole role, final int mapq)
        {
            SAMRecord record = record(readName, role);
            assertNotNull("no emitted record for " + readName + "/" + role, record);
            assertEquals(readName + "/" + role + " mapq", mapq, record.getMappingQuality());
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
