package com.hartwig.hmftools.purple.config;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface DBConfig {

    String DB_ENABLED = "db_enabled";
    String DB_USER = "db_user";
    String DB_PASS = "db_pass";
    String DB_URL = "db_url";

    static void addOptions(@NotNull Options options) {
        options.addOption(DB_ENABLED, false, "Persist data to DB.");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url in form: mysql://host:port/database");
    }

    boolean enabled();

    @NotNull
    String user();

    @NotNull
    String password();

    @NotNull
    String url();

    @NotNull
    static DBConfig createConfig(@NotNull final CommandLine cmd) {
        boolean enabled = cmd.hasOption(DB_ENABLED);
        return ImmutableDBConfig.builder()
                .enabled(enabled)
                .user(enabled ? cmd.getOptionValue(DB_USER) : "")
                .password(enabled ? cmd.getOptionValue(DB_PASS) : "")
                .url("jdbc:" + (enabled ? cmd.getOptionValue(DB_URL) : ""))
                .build();
    }
}
