package Utils;

import static java.util.stream.Collectors.toList;

import static Utils.SvTestRoutines.initialiseSV;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;

import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFile;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.jetbrains.annotations.NotNull;

public class SampleDataLoader
{
    public static List<StructuralVariantData> fromResource(@NotNull final String resource)
    {
        final InputStream inputStream = StructuralVariantFile.class.getResourceAsStream("/sample_data/" + resource);
        return new BufferedReader(new InputStreamReader(inputStream)).lines()
                .filter(x -> !x.startsWith("svId")).map(StructuralVariantFile::fromString).collect(toList());
    }

    public static List<SvVarData> createSVs(final List<StructuralVariantData> svDataList)
    {
        List<SvVarData> svList = svDataList.stream().map(x -> new SvVarData(x)).collect(toList());
        svList.forEach(x -> initialiseSV(x));
        return svList;
    }

}
