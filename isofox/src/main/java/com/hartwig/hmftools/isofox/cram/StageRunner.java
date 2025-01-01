package com.hartwig.hmftools.isofox.cram;

import com.google.auth.oauth2.GoogleCredentials;
import com.google.cloud.http.HttpTransportOptions;
import com.google.cloud.storage.Storage;
import com.google.cloud.storage.StorageOptions;
import com.hartwig.computeengine.execution.vm.BashStartupScript;
import com.hartwig.computeengine.execution.vm.RuntimeFiles;
import com.hartwig.computeengine.execution.vm.command.BashCommand;
import com.hartwig.computeengine.execution.vm.command.OutputUploadCommand;
import com.hartwig.computeengine.storage.GoogleStorageLocation;
import com.hartwig.computeengine.storage.RuntimeBucket;
import com.hartwig.computeengine.storage.RuntimeBucketOptions;

import java.io.IOException;
import java.util.List;


public class StageRunner {

    public void run(String inputPath, String outputPath) throws IOException {
        final List<BashCommand> commands = new CramAndValidateCommands(inputPath, outputPath).commands();

        GoogleCredentials credentials = GoogleCredentials.getApplicationDefault();

        StorageOptions.Builder builder = StorageOptions.newBuilder();
        Storage storage = builder.setCredentials(credentials)
                .setProjectId("hmf-pipeline-development")
                .setTransportOptions(HttpTransportOptions.newBuilder()
                        .setConnectTimeout(0)
                        .setReadTimeout(0)
                        .build())
                .build()
                .getService();


//        RunIdentifier runIdentifier = ArgumentUtil.toRunIdentifier(arguments, metadata);
        var runtimeBucketOptions = RuntimeBucketOptions.builder()
                .namespace("pilot-1")
                .region("europe-west4")
//                .labels(labels.asMap())
//                .runIdentifier(runIdentifier)
                .cmek("projects/hmf-pipeline-development/locations/europe-west4/keyRings" + "/hmf-pipeline-development/cryptoKeys/default-test")
                .build();
        RuntimeBucket bucket = RuntimeBucket.from(storage, runtimeBucketOptions);

        BashStartupScript bash = BashStartupScript.of(bucket.name());
        bash.addCommands(commands)
                .addCommand(new OutputUploadCommand(GoogleStorageLocation.of(bucket.name(), outputPath),
                        RuntimeFiles.typical()));
    }


}
