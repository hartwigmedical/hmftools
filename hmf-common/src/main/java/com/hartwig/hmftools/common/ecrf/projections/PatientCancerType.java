package com.hartwig.hmftools.common.ecrf.projections;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import com.google.common.collect.Lists;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVPrinter;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class PatientCancerType {
    private enum Header {
        patientIdentifier,
        primaryTumorLocation,
        cancerSubtype
    }

    @NotNull
    public abstract String patientIdentifier();

    @NotNull
    public abstract String primaryTumorLocation();

    @NotNull
    public abstract String cancerSubtype();

    @NotNull
    private List<String> csvRecord() {
        return Lists.newArrayList(patientIdentifier(), primaryTumorLocation(), cancerSubtype());
    }

    public static void writeRecords(@NotNull final String outputPath, @NotNull final List<PatientCancerType> patientCancerTypes)
            throws IOException {
        final CSVFormat format = CSVFormat.DEFAULT.withNullString("").withHeader(Header.class);
        final CSVPrinter printer = new CSVPrinter(new FileWriter(outputPath), format);
        printer.printRecords(patientCancerTypes.stream().map(PatientCancerType::csvRecord).collect(Collectors.toList()));
        printer.close();
    }

    @NotNull
    public static List<PatientCancerType> readRecords(@NotNull final String filePath) throws IOException {
        final CSVParser parser = CSVParser.parse(new File(filePath),
                Charset.defaultCharset(),
                CSVFormat.DEFAULT.withHeader(Header.class).withSkipHeaderRecord());
        return StreamSupport.stream(parser.spliterator(), false)
                .map(record -> ImmutablePatientCancerType.of(record.get(Header.patientIdentifier),
                        record.get(Header.primaryTumorLocation),
                        record.get(Header.cancerSubtype)))
                .collect(Collectors.toList());
    }
}
