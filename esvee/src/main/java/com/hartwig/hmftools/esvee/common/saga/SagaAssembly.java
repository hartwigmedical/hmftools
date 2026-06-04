package com.hartwig.hmftools.esvee.common.saga;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;

import java.util.List;

public record SagaAssembly(
        String fastaLabel,
        SagaVariant variant,
        // Junctions occur just before each of these indices in the assembly sequence.
        // Usually length 2. Can be length 1 for DELs where the junction is a single position.
        List<Integer> junctionOffsets,
        int assemblyLength
)
{
    public SagaAssembly
    {
        // E.g.
        // sequence = RRRJJJRRR
        // junctionOffset[0] = 3
        // junctionOffset[1] = 6
        if(!junctionOffsets.stream().allMatch(offset -> offset >= 1 && offset < assemblyLength))
        {
            throw new IllegalArgumentException("Junction offsets out of bounds");
        }
    }

    public String variantId()
    {
        return variant.id();
    }

    public static SagaAssembly fromFastaLabel(final String fastaLabel, int assemblyLength)
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
        SagaBreakend breakend1 = SagaBreakend.fromString(parts[1]);
        SagaBreakend breakend2 = SagaBreakend.fromString(parts[2]);
        int junctionOffset1 = Integer.parseInt(parts[3]);
        Integer junctionOffset2 = parts.length >= 5 ? Integer.parseInt(parts[4]) : null;
        if(junctionOffset2 != null && junctionOffset1 >= junctionOffset2)
        {
            throw new IllegalArgumentException("Invalid junction offsets");
        }
        SagaVariant variant = new SagaVariant(id, breakend1, breakend2);
        List<Integer> junctionOffsets = junctionOffset2 == null ? List.of(junctionOffset1) : List.of(junctionOffset1, junctionOffset2);
        return new SagaAssembly(fastaLabel, variant, junctionOffsets, assemblyLength);
    }
}
