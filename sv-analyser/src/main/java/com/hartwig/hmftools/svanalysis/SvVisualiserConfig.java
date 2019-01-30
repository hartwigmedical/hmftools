package com.hartwig.hmftools.svanalysis;

import static java.util.stream.Collectors.toList;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.svanalysis.visualisation.CopyNumberAlteration;
import com.hartwig.hmftools.svanalysis.visualisation.CopyNumberAlterations;
import com.hartwig.hmftools.svanalysis.visualisation.ImmutableCopyNumberAlteration;
import com.hartwig.hmftools.svanalysis.visualisation.Link;
import com.hartwig.hmftools.svanalysis.visualisation.Links;
import com.hartwig.hmftools.svanalysis.visualisation.Segment;
import com.hartwig.hmftools.svanalysis.visualisation.Segments;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface SvVisualiserConfig {

    Logger LOGGER = LogManager.getLogger(SvVisualiserConfig.class);

    String PLOT_OUT = "plot_out";
    String DATA_OUT = "data_out";
    String SAMPLE = "sample";
    String SEGMENT = "segment";
    String LINK = "link";
    String CIRCOS = "circos";
    String THREADS = "threads";
    String DEBUG = "debug";
    String SINGLE_CLUSTER = "clusterId";
    String SINGLE_CHROMOSOME = "chromosome";
    String DB_USER = "db_user";
    String DB_PASS = "db_pass";
    String DB_URL = "db_url";
    String CNA = "cna";

    @NotNull
    String sample();

    @NotNull
    List<Segment> segments();

    @NotNull
    List<Link> links();

    @NotNull
    List<CopyNumberAlteration> copyNumberAlterations();

    @NotNull
    String outputConfPath();

    @NotNull
    String outputPlotPath();

    @NotNull
    String circosBin();

    int threads();

    boolean debug();

    @Nullable
    Integer singleCluster();

    @Nullable
    String singleChromosome();

    @NotNull
    static Options createOptions() {
        final Options options = new Options();
        options.addOption(PLOT_OUT, true, "Plot output directory");
        options.addOption(DATA_OUT, true, "Data output directory");
        options.addOption(SAMPLE, true, "Sample name");
        options.addOption(SEGMENT, true, "Path to track file");
        options.addOption(LINK, true, "Path to link file");
        options.addOption(CIRCOS, true, "Path to circos binary");
        options.addOption(THREADS, true, "Number of threads to use");
        options.addOption(DEBUG, false, "Enabled debug mode");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url in form: mysql://host:port/database");
        options.addOption(SINGLE_CLUSTER, true, "Only generate image for single cluster");
        options.addOption(SINGLE_CHROMOSOME, true, "Only generate image for singe chromosome");
        options.addOption(CNA, true, "Location of copy number alterations (optional alternative to db)");

        return options;
    }

    @NotNull
    static SvVisualiserConfig createConfig(@NotNull final CommandLine cmd) throws ParseException, IOException, SQLException {
        final StringJoiner missingJoiner = new StringJoiner(", ");
        final String linkPath = parameter(cmd, LINK, missingJoiner);
        final String trackPath = parameter(cmd, SEGMENT, missingJoiner);
        final String sample = parameter(cmd, SAMPLE, missingJoiner);
        final String plotOutputDir = parameter(cmd, PLOT_OUT, missingJoiner);
        final String dataOutputDir = parameter(cmd, DATA_OUT, missingJoiner);
        final String circos = parameter(cmd, CIRCOS, missingJoiner);
        final String dbUser;
        final String dbPassword;
        final String dbUrl;
        if (!cmd.hasOption(CNA)) {
            dbUser = parameter(cmd, DB_USER, missingJoiner);
            dbPassword = parameter(cmd, DB_PASS, missingJoiner);
            dbUrl = parameter(cmd, DB_URL, missingJoiner);
        } else {
            dbUser = "";
            dbPassword = "";
            dbUrl = "";
        }

        final String missing = missingJoiner.toString();

        if (!missing.isEmpty()) {
            throw new ParseException("Missing the following parameters: " + missing);
        }

        final List<Segment> segments =
                Segments.readTracksFromFile(trackPath).stream().filter(x -> x.sampleId().equals(sample)).collect(toList());
        final List<Link> links = Links.readLinks(linkPath).stream().filter(x -> x.sampleId().equals(sample)).collect(toList());

        if (segments.isEmpty() && links.isEmpty()) {
            LOGGER.warn("No structural variants found for sample {}", sample);
        }

        final List<CopyNumberAlteration> cna;
        if (cmd.hasOption(CNA)) {
            LOGGER.info("Reading copy numbers from {}", cmd.getOptionValue(CNA));
            cna = CopyNumberAlterations.read(cmd.getOptionValue(CNA));
        } else {
            LOGGER.info("Loading copy numbers from database");
            cna = sampleCopyNumberAlterations(sample, dbUser, dbPassword, "jdbc:" + dbUrl);
        }

        if (cna.isEmpty()) {
            LOGGER.warn("No copy number alterations found for sample {}", sample);
        }

        return ImmutableSvVisualiserConfig.builder()
                .outputConfPath(dataOutputDir)
                .outputPlotPath(plotOutputDir)
                .segments(segments)
                .links(links)
                .sample(sample)
                .copyNumberAlterations(cna)
                .circosBin(circos)
                .threads(Integer.valueOf(cmd.getOptionValue(THREADS, "1")))
                .debug(cmd.hasOption(DEBUG))
                .singleCluster(cmd.hasOption(SINGLE_CLUSTER) ? Integer.valueOf(cmd.getOptionValue(SINGLE_CLUSTER)) : null)
                .singleChromosome(cmd.hasOption(SINGLE_CHROMOSOME) ? cmd.getOptionValue(SINGLE_CHROMOSOME) : null)
                .build();
    }

    @NotNull
    static String parameter(@NotNull final CommandLine cmd, @NotNull final String parameter, @NotNull final StringJoiner missing) {
        final String value = cmd.getOptionValue(parameter);
        if (value == null) {
            missing.add(parameter);
            return "";
        }
        return value;
    }

    @NotNull
    static List<CopyNumberAlteration> sampleCopyNumberAlterations(@NotNull final String sample, @NotNull final String userName,
            @NotNull final String password, @NotNull final String url) throws SQLException {
        final DatabaseAccess dbAccess = new DatabaseAccess(userName, password, url);
        return dbAccess.readCopynumbers(sample)
                .stream()
                .map(x -> ImmutableCopyNumberAlteration.builder()
                        .from(x)
                        .baf(x.averageActualBAF())
                        .copyNumber(x.averageTumorCopyNumber())
                        .build())
                .collect(toList());
    }

}
