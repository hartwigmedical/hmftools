
plugins {
    id("com.hartwig.java-conventions")
}

description = "HMF Tools - PURPLE"

dependencies {
    implementation(project(":hmf-common"))
    implementation(project(":patient-db"))
    annotationProcessor(libs.immutables.value)
    compileOnly(libs.immutables.value)
    testImplementation(libs.mockito.core)
    testImplementation(libs.junit)
}

shadowJar {
    mainClass.set("com.hartwig.hmftools.purple.PurpleApplication")
}
