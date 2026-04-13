package com.hartwig.hmftools.esvee.common;

import static java.util.function.UnaryOperator.identity;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BasePosition;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class SagaResource
{
    private final String mFastaPath;
    private final Map<String, AssemblyMetadata> mAssembliesById;
    private final Map<String, List<IndexedBreakend>> mSearchableBreakends;
    private SAMSequenceDictionary mSamDict;

    public SagaResource(final String fastaPath)
    {
        mFastaPath = fastaPath;
        List<AssemblyMetadata> assemblies = loadAssemblies();
        mAssembliesById = assemblies.stream()
                .collect(Collectors.toMap(AssemblyMetadata::id, identity()));
        mSearchableBreakends = assemblies.stream()
                .flatMap(assembly ->
                        assembly.variant.breakends().map(breakend ->
                                new IndexedBreakend(breakend, assembly.variant.id())))
                .collect(Collectors.groupingBy(IndexedBreakend::chromosome));
    }

    public String bwaIndexImagePath()
    {
        return mFastaPath + ".img";
    }

    public SAMSequenceDictionary samDict()
    {
        return mSamDict;
    }

    public AssemblyMetadata getAssemblyById(final String variantId)
    {
        AssemblyMetadata assembly = mAssembliesById.get(variantId);
        if(assembly == null)
        {
            throw new IllegalArgumentException("No SAGA variant with ID:" + variantId);
        }
        else
        {
            return assembly;
        }
    }

    public Variant getVariantById(final String variantId)
    {
        return getAssemblyById(variantId).variant();
    }

    public AssemblyMetadata getAssemblyByFastaLabel(final String fastaLabel)
    {
        String variantId = fastaLabel.split("\\|", 2)[0];
        return getAssemblyById(variantId);
    }

    public Map<String, List<IndexedBreakend>> searchableBreakends()
    {
        return mSearchableBreakends;
    }

    public record Breakend(
            BasePosition position,
            Orientation orientation
    )
    {
        public static Breakend fromString(final String string)
        {
            String[] parts = string.split(":");
            if(parts.length != 3)
            {
                throw new IllegalArgumentException("Invalid breakend string: " + string);
            }
            String chromosome = parts[0];
            int position = Integer.parseInt(parts[1]);
            Orientation orientation = Orientation.fromByteStr(parts[2]);
            return new Breakend(new BasePosition(chromosome, position), orientation);
        }
    }

    public record Variant(
            String id,
            Breakend breakend1,
            Breakend breakend2
    )
    {
        public Stream<Breakend> breakends()
        {
            return Stream.of(breakend1, breakend2);
        }
    }

    public record AssemblyMetadata(
            String fastaLabel,
            Variant variant,
            // Usually length 2. Can be length 1 for DELs where the junction is a single position.
            List<Integer> junctionOffsets
    )
    {
        public String id()
        {
            return variant.id();
        }

        public static AssemblyMetadata fromFastaLabel(final String fastaLabel)
        {
            // E.g.:
            // SvimAsm00000237|chr1:181626:1|chr1:181627:-1|150|285
            // SvimAsm00000238|chr1:368909:1|chr1:369380:-1|150|

            String[] parts = fastaLabel.split("\\|");
            if(!(parts.length == 4 || parts.length == 5))
            {
                SV_LOGGER.error("Expected 4 or 5 parts but got {}", parts.length);
                throw new IllegalArgumentException("Invalid fasta label: " + fastaLabel);
            }
            String id = parts[0];
            Breakend breakend1 = Breakend.fromString(parts[1]);
            Breakend breakend2 = Breakend.fromString(parts[2]);
            int junctionOffset1 = Integer.parseInt(parts[3]);
            Integer junctionOffset2 = parts.length >= 5 ? Integer.parseInt(parts[4]) : null;
            if(junctionOffset2 != null && junctionOffset1 >= junctionOffset2)
            {
                throw new IllegalArgumentException("Invalid junction offsets");
            }
            Variant variant = new Variant(id, breakend1, breakend2);
            List<Integer> junctionOffsets = junctionOffset2 == null ? List.of(junctionOffset1) : List.of(junctionOffset1, junctionOffset2);
            return new AssemblyMetadata(fastaLabel, variant, junctionOffsets);
        }
    }

    public record IndexedBreakend(
            Breakend breakend,
            String variantId
    )
    {
        public BasePosition position()
        {
            return breakend.position();
        }

        public String chromosome()
        {
            return position().Chromosome;
        }
    }

    private List<AssemblyMetadata> loadAssemblies()
    {
        SV_LOGGER.debug("Loading SAGA resource FASTA");
        List<AssemblyMetadata> assemblies;
        try(IndexedFastaSequenceFile fasta = loadFasta())
        {
            mSamDict = fasta.getSequenceDictionary();
            assemblies = mSamDict.getSequences().stream()
                    .map(seq -> AssemblyMetadata.fromFastaLabel(seq.getSequenceName()))
                    .toList();
        }
        catch(IOException e)
        {
            // Failed to load/read/close fasta file.
            throw new RuntimeException("Failed to load SAGA resource", e);
        }

        // Variants are keyed by ID later, so ensure the IDs are unique.
        long uniqueVariantIdCount = assemblies.stream().map(assembly -> assembly.variant().id()).distinct().count();
        if(uniqueVariantIdCount != assemblies.size())
        {
            throw new RuntimeException("Duplicate variant IDs in SAGA resource");
        }

        SV_LOGGER.debug("Loaded {} SAGA variant assemblies", assemblies.size());

        return assemblies;
    }

    private IndexedFastaSequenceFile loadFasta() throws IOException
    {
        return new IndexedFastaSequenceFile(new File(mFastaPath));
    }
}
