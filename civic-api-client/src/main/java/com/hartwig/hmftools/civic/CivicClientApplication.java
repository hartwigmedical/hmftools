package com.hartwig.hmftools.civic;

import com.google.gson.Gson;
import com.hartwig.hmftools.civic.api.CivicApi;
import com.hartwig.hmftools.civic.data.CivicApiDataGson;
import com.hartwig.hmftools.civic.data.CivicGene;
import com.hartwig.hmftools.civic.data.CivicVariant;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import okhttp3.OkHttpClient;
import retrofit2.Retrofit;
import retrofit2.adapter.rxjava2.RxJava2CallAdapterFactory;
import retrofit2.converter.gson.GsonConverterFactory;

public class CivicClientApplication {

    private static final Logger LOGGER = LogManager.getLogger(CivicClientApplication.class);

    private static final String CIVIC_API_ENDPOINT = "https://civic.genome.wustl.edu/api/";

    public static void main(final String... args) throws Exception {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        OkHttpClient client = new OkHttpClient.Builder().build();

        final Gson gson = CivicApiDataGson.buildGson();

        final Retrofit retrofit = new Retrofit.Builder().baseUrl(CIVIC_API_ENDPOINT)
                .addConverterFactory(GsonConverterFactory.create(gson))
                .addCallAdapterFactory(RxJava2CallAdapterFactory.create())
                .client(client)
                .build();

        final CivicApi civic = retrofit.create(CivicApi.class);
        civic.getGene(1956)
                .flatMapIterable(CivicGene::variants)
                .filter(variantKey -> variantKey.acceptedEvidenceItems() > 0)
                .flatMap(variantKey -> civic.getVariant(variantKey.id()))
                .flatMapIterable(CivicVariant::evidenceItems)
                .filter(evidenceItem -> !evidenceItem.drugs().isEmpty())
                .filter(evidenceItem -> evidenceItem.level() <= 'C')
                .subscribe(LOGGER::info, Throwable::printStackTrace, () -> LOGGER.info("completed"));
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
