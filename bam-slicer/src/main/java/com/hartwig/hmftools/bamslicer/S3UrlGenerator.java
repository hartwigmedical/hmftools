package com.hartwig.hmftools.bamslicer;

import java.io.IOException;
import java.net.URL;
import java.util.Date;

import com.amazonaws.AmazonClientException;
import com.amazonaws.AmazonServiceException;
import com.amazonaws.HttpMethod;
import com.amazonaws.auth.profile.ProfileCredentialsProvider;
import com.amazonaws.client.builder.AwsClientBuilder;
import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.AmazonS3ClientBuilder;
import com.amazonaws.services.s3.model.GeneratePresignedUrlRequest;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
abstract class S3UrlGenerator {
    private static final Logger LOGGER = LogManager.getLogger(S3UrlGenerator.class);

    abstract String endpoint_url();

    abstract String profile();

    @Value.Lazy
    AmazonS3 s3Client() {
        return AmazonS3ClientBuilder.standard()
                .withEndpointConfiguration(new AwsClientBuilder.EndpointConfiguration(endpoint_url(), null))
                .withCredentials(new ProfileCredentialsProvider(profile()))
                .build();
    }

    @NotNull
    URL generateUrl(@NotNull final String bucketName, @NotNull final String objectKey, final int expirationHours) throws IOException {
        try {
            LOGGER.info("Generating pre-signed URL for bucket: {}\tobject: {}\texpirationTime: {} hours",
                    bucketName,
                    objectKey,
                    expirationHours);
            final long millisNow = System.currentTimeMillis();
            final long expirationMillis = millisNow + 1000 * 60 * 60 * expirationHours;
            final Date expiration = new java.util.Date(expirationMillis);
            final GeneratePresignedUrlRequest generatePresignedUrlRequest = new GeneratePresignedUrlRequest(bucketName, objectKey);
            generatePresignedUrlRequest.setMethod(HttpMethod.GET);
            generatePresignedUrlRequest.setExpiration(expiration);

            return s3Client().generatePresignedUrl(generatePresignedUrlRequest);
        } catch (AmazonServiceException exception) {
            LOGGER.error("Error Message: {}", exception.getMessage());
            LOGGER.error("HTTP  Code: {}", exception.getStatusCode());
            LOGGER.error("Error Code: {}", exception.getErrorCode());
            LOGGER.error("Error Type:   {}", exception.getErrorType());
            LOGGER.error("Request ID:   {}", exception.getRequestId());
        } catch (AmazonClientException ace) {
            LOGGER.error("The client encountered an internal error while trying to communicate with S3");
            LOGGER.error("Error Message: {}", ace.getMessage());
        }
        throw new IOException("Failed to create sbp URL.");
    }
}