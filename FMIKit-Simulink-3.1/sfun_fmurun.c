/*****************************************************************
 *  Copyright (c) Dassault Systemes. All rights reserved.        *
 *  This file is part of FMIKit. See LICENSE.txt in the project  *
 *  root for license information.                                *
 *****************************************************************/

#define FMI_MAX_MESSAGE_LENGTH 4096

#define INTERNET_MAX_URL_LENGTH 2083

#ifdef GRTFMI
extern const char *FMU_RESOURCES_DIR;
#endif

#define FMI_STATIC static

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>

#ifdef _WIN32
#include <shlwapi.h>
#pragma comment(lib, "shlwapi.lib")
#else
#include <stdarg.h>
#include <dlfcn.h>
#endif

#ifdef _MSC_VER
#define strdup _strdup
#endif

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <Windows.h>
#endif

#include <stdbool.h>
#include <stddef.h>

#ifndef FMI_MAX_MESSAGE_LENGTH
#define FMI_MAX_MESSAGE_LENGTH 4096
#endif

#ifndef FMI_STATIC
#define FMI_STATIC
#endif

typedef enum {
    FMIOK,
    FMIWarning,
    FMIDiscard,
    FMIError,
    FMIFatal,
    FMIPending
} FMIStatus;

typedef enum {

    // FMI 3.0 variable types
    FMIFloat32Type,
    FMIDiscreteFloat32Type,
    FMIFloat64Type,
    FMIDiscreteFloat64Type,
    FMIInt8Type,
    FMIUInt8Type,
    FMIInt16Type,
    FMIUInt16Type,
    FMIInt32Type,
    FMIUInt32Type,
    FMIInt64Type,
    FMIUInt64Type,
    FMIBooleanType,
    FMIStringType,
    FMIBinaryType,
    FMIClockType,

    // non-variable types used as parameters
    FMIValueReferenceType,
    FMISizeTType,

    // aliases for FMI 1.0 and 2.0 variable types
    FMIRealType         = FMIFloat64Type,
    FMIDiscreteRealType = FMIDiscreteFloat64Type,
    FMIIntegerType      = FMIInt32Type,

} FMIVariableType;

typedef enum {
    FMIVersion1 = 1,
    FMIVersion2 = 2,
    FMIVersion3 = 3
} FMIVersion;

typedef enum {
    FMIModelExchange,
    FMICoSimulation,
    FMIScheduledExecution
} FMIInterfaceType;

typedef enum {
    FMI2StartAndEndState        = 1 << 0,
    FMI2InstantiatedState       = 1 << 1,
    FMI2InitializationModeState = 1 << 2,

    // model exchange states
    FMI2EventModeState          = 1 << 3,
    FMI2ContinuousTimeModeState = 1 << 4,

    // co-simulation states
    FMI2StepCompleteState       = 1 << 5,
    FMI2StepInProgressState     = 1 << 6,
    FMI2StepFailedState         = 1 << 7,
    FMI2StepCanceledState       = 1 << 8,

    FMI2TerminatedState         = 1 << 9,
    FMI2ErrorState              = 1 << 10,
    FMI2FatalState              = 1 << 11,
} FMI2State;

typedef unsigned int FMIValueReference;

typedef struct FMIInstance_ FMIInstance;

typedef struct FMI1Functions_ FMI1Functions;

typedef struct FMI2Functions_ FMI2Functions;

typedef struct FMI3Functions_ FMI3Functions;

typedef void FMILogFunctionCall(FMIInstance *instance, FMIStatus status, const char *message);

typedef void FMILogMessage(FMIInstance *instance, FMIStatus status, const char *category, const char *message);

struct FMIInstance_ {

    FMI1Functions *fmi1Functions;
    FMI2Functions *fmi2Functions;
    FMI3Functions *fmi3Functions;

#ifdef _WIN32
    HMODULE libraryHandle;
#else
    void *libraryHandle;
#endif

    void *userData;

    FMILogMessage      *logMessage;
    FMILogFunctionCall *logFunctionCall;

    double time;

    char* logMessageBuffer;
    size_t logMessageBufferSize;
    size_t logMessageBufferPosition;

    void *component;

    const char *name;

    bool logFMICalls;

    FMI2State state;

    FMIStatus status;

    FMIVersion fmiVersion;

    FMIInterfaceType interfaceType;

};

FMI_STATIC FMIInstance* FMICreateInstance(const char* instanceName, const char* libraryPath, FMILogMessage* logMessage, FMILogFunctionCall* logFunctionCall);

FMI_STATIC void FMIFreeInstance(FMIInstance *instance);

FMI_STATIC void FMIClearLogMessageBuffer(FMIInstance* instance);

FMI_STATIC void FMIAppendToLogMessageBuffer(FMIInstance* instance, const char* format, ...);

FMI_STATIC void FMIAppendArrayToLogMessageBuffer(FMIInstance* instance, const void* values, size_t nValues, const size_t sizes[], FMIVariableType variableType);

FMI_STATIC FMIStatus FMIURIToPath(const char *uri, char *path, const size_t pathLength);

FMI_STATIC FMIStatus FMIPathToURI(const char *path, char *uri, const size_t uriLength);

FMI_STATIC FMIStatus FMIPlatformBinaryPath(const char *unzipdir, const char *modelIdentifier, FMIVersion fmiVersion, char *platformBinaryPath, size_t size);

#ifdef __cplusplus
}  /* end of extern "C" { */
#endif


FMIInstance *FMICreateInstance(const char *instanceName, const char *libraryPath, FMILogMessage *logMessage, FMILogFunctionCall *logFunctionCall) {

# ifdef _WIN32
    TCHAR Buffer[1024];
    GetCurrentDirectory(1024, Buffer);

    WCHAR dllDirectory[MAX_PATH];

    // convert path to unicode
    mbstowcs(dllDirectory, libraryPath, MAX_PATH);

    // remove the file name
    PathRemoveFileSpecW(dllDirectory);

    // add the binaries directory temporarily to the DLL path to allow discovery of dependencies
    DLL_DIRECTORY_COOKIE dllDirectoryCookie = AddDllDirectory(dllDirectory);

    // TODO: log getLastSystemError()

    HMODULE libraryHandle = LoadLibraryExA(libraryPath, NULL, LOAD_LIBRARY_SEARCH_DEFAULT_DIRS);

    // remove the binaries directory from the DLL path
    if (dllDirectoryCookie) {
        RemoveDllDirectory(dllDirectoryCookie);
    }

    // TODO: log error

# else
    void *libraryHandle = dlopen(libraryPath, RTLD_LAZY);
# endif

    if (!libraryHandle) {
        return NULL;
    }

    FMIInstance* instance = (FMIInstance*)calloc(1, sizeof(FMIInstance));

    instance->libraryHandle = libraryHandle;

    instance->logMessage = logMessage;
    instance->logFunctionCall = logFunctionCall;

    instance->logMessageBufferSize = 1024;
    instance->logMessageBuffer = (char*)calloc(instance->logMessageBufferSize, sizeof(char));
    instance->logMessageBufferPosition = 0;

    instance->name = strdup(instanceName);

    instance->status = FMIOK;

    return instance;
}

void FMIFreeInstance(FMIInstance *instance) {

    if (!instance) {
        return;
    }

    // unload the shared library
    if (instance->libraryHandle) {
# ifdef _WIN32
        FreeLibrary(instance->libraryHandle);
# else
        dlclose(instance->libraryHandle);
# endif
        instance->libraryHandle = NULL;
    }

    free(instance->logMessageBuffer);

    free((void*)instance->name);

    free(instance->fmi1Functions);
    free(instance->fmi2Functions);
    free(instance->fmi3Functions);

    free(instance);
}

void FMIClearLogMessageBuffer(FMIInstance* instance) {
    instance->logMessageBufferPosition = 0;
    snprintf(instance->logMessageBuffer, instance->logMessageBufferSize, "");
}

void FMIAppendToLogMessageBuffer(FMIInstance* instance, const char* format, ...) {

    va_list args;

    va_start(args, format);
    const int length = vsnprintf(&instance->logMessageBuffer[instance->logMessageBufferPosition], instance->logMessageBufferSize - instance->logMessageBufferPosition, format, args);
    va_end(args);

    if (length < instance->logMessageBufferSize - instance->logMessageBufferPosition) {

        instance->logMessageBufferPosition += length;

    } else {

        while (instance->logMessageBufferSize < instance->logMessageBufferPosition + length) {
            instance->logMessageBufferSize *= 2;
        }

        instance->logMessageBuffer = realloc(instance->logMessageBuffer, instance->logMessageBufferSize);

        va_start(args, format);
        instance->logMessageBufferPosition += vsnprintf(&instance->logMessageBuffer[instance->logMessageBufferPosition], instance->logMessageBufferSize - instance->logMessageBufferPosition, format, args);
        va_end(args);

    }
}

void FMIAppendArrayToLogMessageBuffer(FMIInstance* instance, const void* values, size_t nValues, const size_t sizes[], FMIVariableType variableType) {

    for (size_t i = 0; i < nValues;) {

        // pointer to the last byte (terminator)
        char* s = &instance->logMessageBuffer[instance->logMessageBufferPosition];

        // remaining bytes in the buffer
        size_t n = instance->logMessageBufferSize - instance->logMessageBufferPosition;

        int length;

        switch (variableType) {
        case FMIFloat32Type:
        case FMIDiscreteFloat32Type:
            length = snprintf(s, n, "%.7g", ((float*)values)[i]);
            break;
        case FMIFloat64Type:
        case FMIDiscreteFloat64Type:
            length = snprintf(s, n, "%.16g", ((double*)values)[i]);
            break;
        case FMIInt8Type:
            length = snprintf(s, n, "%" PRId8, ((int8_t *)values)[i]);
            break;
        case FMIUInt8Type:
            length = snprintf(s, n, "%" PRIu8, ((uint8_t *)values)[i]);
            break;
        case FMIInt16Type:
            length = snprintf(s, n, "%" PRId16, ((int16_t *)values)[i]);
            break;
        case FMIUInt16Type:
            length = snprintf(s, n, "%" PRIu16, ((uint16_t *)values)[i]);
            break;
        case FMIInt32Type:
            length = snprintf(s, n, "%" PRId32, ((int32_t *)values)[i]);
            break;
        case FMIUInt32Type:
            length = snprintf(s, n, "%" PRIu32, ((uint32_t*)values)[i]);
            break;
        case FMIInt64Type:
            length = snprintf(s, n, "%" PRId64, ((int64_t *)values)[i]);
            break;
        case FMIUInt64Type:
            length = snprintf(s, n, "%" PRIu64, ((uint64_t *)values)[i]);
            break;
        case FMIBooleanType:
            switch (instance->fmiVersion) {
                case FMIVersion1:
                    length = snprintf(s, n, "%d", ((char*)values)[i]);
                    break;
                case FMIVersion2:
                    length = snprintf(s, n, "%d", ((int*)values)[i]);
                    break;
                case FMIVersion3:
                    length = snprintf(s, n, "%d", ((bool*)values)[i]);
                    break;
            }
            break;
        case FMIStringType:
            length = snprintf(s, n, "\"%s\"", ((const char**)values)[i]);
            break;
        case FMIBinaryType: {
            const size_t size = sizes[i];
            const unsigned char* v = ((const unsigned char**)values)[i];
            length = 2 + size * 2;
            if (length < n) {
                snprintf(s, n, "0x");
                s += 2;
                n -= 2;
                for (size_t j = 0; j < size; j++) {
                    snprintf(s, n, "%02hhx", v[j]);
                    s += 2;
                    n -= 2;
                }
            }
            break;
        }
        case FMIClockType:
            length = snprintf(s, n, "%d", ((bool *)values)[i]);
            break;
        case FMIValueReferenceType:
            length = snprintf(s, n, "%u", ((FMIValueReference*)values)[i]);
            break;
        case FMISizeTType:
            length = snprintf(s, n, "%zu", ((size_t*)values)[i]);
            break;
        default:
            i++;
            continue;
        }

        if (length + sizeof(", ") < instance->logMessageBufferSize - instance->logMessageBufferPosition) {

            instance->logMessageBufferPosition += length;
            i++;

        } else {

            while (instance->logMessageBufferSize < instance->logMessageBufferPosition + length) {
                instance->logMessageBufferSize *= 2;
            }

            instance->logMessageBuffer = realloc(instance->logMessageBuffer, instance->logMessageBufferSize);

            continue;
        }

        if (i < nValues) {
            s = &instance->logMessageBuffer[instance->logMessageBufferPosition];
            n = instance->logMessageBufferSize - instance->logMessageBufferPosition;
            instance->logMessageBufferPosition += snprintf(s, n, ", ");
        }
    }
}

FMIStatus FMIURIToPath(const char *uri, char *path, const size_t pathLength) {

#ifdef _WIN32
    DWORD pcchPath = (DWORD)pathLength;

    if (PathCreateFromUrlA(uri, path, &pcchPath, 0) != S_OK) {
        return FMIError;
    }
#else
    const char *scheme1 = "file:///";
    const char *scheme2 = "file:/";

    strncpy(path, uri, pathLength);

    if (strncmp(uri, scheme1, strlen(scheme1)) == 0) {
        strncpy(path, &uri[strlen(scheme1)] - 1, pathLength);
    } else if (strncmp(uri, scheme2, strlen(scheme2)) == 0) {
        strncpy(path, &uri[strlen(scheme2) - 1], pathLength);
    } else {
        return FMIError;
    }
#endif

#ifdef _WIN32
    const char* sep = "\\";
#else
    const char* sep = "/";
#endif

    if (path[strlen(path) - 1] != sep[0]) {
        strncat(path, sep, pathLength);
    }

    return FMIOK;
}

FMIStatus FMIPathToURI(const char *path, char *uri, const size_t uriLength) {

#ifdef _WIN32
    DWORD pcchUri = (DWORD)uriLength;

    if (UrlCreateFromPathA(path, uri, &pcchUri, 0) != S_OK) {
        return FMIError;
    }
#else
    const size_t pathLen = strlen(path);

    if (uriLength < strlen(path) + 8) {
        return FMIError;
    }

    strcpy(uri, "file://");

    size_t p = 7;

    // percent encode special characters
    for (size_t i = 0; i < pathLen; i++) {

        if (uriLength < p + 4) {
            return FMIError;
        }

        const char c = path[i];

        bool encode = true;

        const char* legal = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-._~:/?#[]@!$&'()*+,;=";

        for (size_t j = 0; j < 84; j++) {
            if (c == legal[j]) {
                encode = false;
                break;
            }
        }

        if (encode) {
            sprintf(&uri[p], "%%%2X", c);
            p += 3;
        } else {
            uri[p++] = c;
        }

    }

    uri[p] = '\0';
#endif

    return FMIOK;
}

FMIStatus FMIPlatformBinaryPath(const char *unzipdir, const char *modelIdentifier, FMIVersion fmiVersion, char *platformBinaryPath, size_t size) {

#if defined(_WIN32)
    const char *platform = "win";
    const char *system   = "windows";
    const char sep       = '\\';
    const char *ext      = ".dll";
#elif defined(__APPLE__)
    const char *platform = "darwin";
    const char *system   = "darwin";
    const char sep       = '/';
    const char *ext      = ".dylib";
#else
    const char *platform = "linux";
    const char *system   = "linux";
    const char sep       = '/';
    const char *ext      = ".so";
#endif

#if defined(_WIN64) || defined(__x86_64__)
    const char *bits = "64";
    const char *arch = "x86_64";
#else
    const char *bits = "32";
    const char *arch = "x86";
#endif
    const char* bin = "binaries";
    char optSep[2] = "";
    int rc;

    if (unzipdir[strlen(unzipdir) - 1] != sep) {
        optSep[0] = sep;
    }
    if (fmiVersion == 3) {
        rc = snprintf(platformBinaryPath, size, "%s%s%s%c%s-%s%c%s%s", unzipdir, optSep, bin, sep, arch, system, sep, modelIdentifier, ext);
    }
    else {
        rc = snprintf(platformBinaryPath, size, "%s%s%s%c%s%s%c%s%s", unzipdir, optSep, bin, sep, platform, bits, sep, modelIdentifier, ext);
    }

    if (rc >= size) {
        return FMIError;
    }

    return FMIOK;
}
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#ifndef _WIN32
#include <dlfcn.h>
#endif

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#ifndef FMI1FUNCTIONS_H
#define FMI1FUNCTIONS_H

/**************************************************************
 *  Copyright (c) Modelica Association Project "FMI".         *
 *  All rights reserved.                                      *
 *  This file is part of the Reference FMUs. See LICENSE.txt  *
 *  in the project root for license information.              *
 **************************************************************/

/* -------------------------------------------------------------------------
 * Combined FMI 1.0 Functions for Model Exchange & Co-Simulation
 * -------------------------------------------------------------------------*/

/* Type definitions */
typedef void*        fmi1Component;
typedef unsigned int fmi1ValueReference;
typedef double       fmi1Real;
typedef int          fmi1Integer;
typedef char         fmi1Boolean;
typedef const char*  fmi1String;

/* Values for fmi1Boolean  */
#define fmi1True  1
#define fmi1False 0

/* Undefined value for fmi1ValueReference (largest unsigned int value) */
#define fmi1UndefinedValueReference (fmi1ValueReference)(-1)


/* make sure all compiler use the same alignment policies for structures */
#if defined _MSC_VER || defined __GNUC__
#pragma pack(push,8)
#endif

typedef enum  {
    fmi1OK,
    fmi1Warning,
    fmi1Discard,
    fmi1Error,
    fmi1Fatal
} fmi1Status;

typedef enum {
    fmi1DoStepStatus,
    fmi1PendingStatus,
    fmi1LastSuccessfulTime
} fmi1StatusKind;

typedef void   (*fmi1CallbackLogger)        (fmi1Component c, fmi1String instanceName, fmi1Status status, fmi1String category, fmi1String message, ...);
typedef void*  (*fmi1CallbackAllocateMemory)(size_t nobj, size_t size);
typedef void   (*fmi1CallbackFreeMemory)    (void* obj);
typedef void   (*fmi1StepFinished)          (fmi1Component c, fmi1Status status);

typedef struct {
    fmi1CallbackLogger         logger;
    fmi1CallbackAllocateMemory allocateMemory;
    fmi1CallbackFreeMemory     freeMemory;
    fmi1StepFinished           stepFinished;
} fmi1CallbackFunctions;

typedef struct {
    fmi1Boolean iterationConverged;
    fmi1Boolean stateValueReferencesChanged;
    fmi1Boolean stateValuesChanged;
    fmi1Boolean terminateSimulation;
    fmi1Boolean upcomingTimeEvent;
    fmi1Real    nextEventTime;
} fmi1EventInfo;

/* reset alignment policy to the one set before reading this file */
#if defined _MSC_VER || defined __GNUC__
#pragma pack(pop)
#endif

/***************************************************
 Common Functions for FMI 1.0
****************************************************/
typedef fmi1Status    fmi1SetRealTYPE         (fmi1Component c, const fmi1ValueReference vr[], size_t nvr, const fmi1Real    value[]);
typedef fmi1Status    fmi1SetIntegerTYPE      (fmi1Component c, const fmi1ValueReference vr[], size_t nvr, const fmi1Integer value[]);
typedef fmi1Status    fmi1SetBooleanTYPE      (fmi1Component c, const fmi1ValueReference vr[], size_t nvr, const fmi1Boolean value[]);
typedef fmi1Status    fmi1SetStringTYPE       (fmi1Component c, const fmi1ValueReference vr[], size_t nvr, const fmi1String  value[]);
typedef fmi1Status    fmi1GetRealTYPE         (fmi1Component c, const fmi1ValueReference vr[], size_t nvr, fmi1Real    value[]);
typedef fmi1Status    fmi1GetIntegerTYPE      (fmi1Component c, const fmi1ValueReference vr[], size_t nvr, fmi1Integer value[]);
typedef fmi1Status    fmi1GetBooleanTYPE      (fmi1Component c, const fmi1ValueReference vr[], size_t nvr, fmi1Boolean value[]);
typedef fmi1Status    fmi1GetStringTYPE       (fmi1Component c, const fmi1ValueReference vr[], size_t nvr, fmi1String  value[]);
typedef fmi1Status    fmi1SetDebugLoggingTYPE (fmi1Component c, fmi1Boolean loggingOn);

/***************************************************
 FMI 1.0 for Model Exchange Functions
****************************************************/
typedef const char*   fmi1GetModelTypesPlatformTYPE      ();
typedef const char*   fmi1GetVersionTYPE                 ();
typedef fmi1Component fmi1InstantiateModelTYPE           (fmi1String instanceName, fmi1String GUID, fmi1CallbackFunctions functions, fmi1Boolean loggingOn);
typedef void          fmi1FreeModelInstanceTYPE          (fmi1Component c);
typedef fmi1Status    fmi1SetTimeTYPE                    (fmi1Component c, fmi1Real time);
typedef fmi1Status    fmi1SetContinuousStatesTYPE        (fmi1Component c, const fmi1Real x[], size_t nx);
typedef fmi1Status    fmi1CompletedIntegratorStepTYPE    (fmi1Component c, fmi1Boolean* callEventUpdate);
typedef fmi1Status    fmi1InitializeTYPE                 (fmi1Component c, fmi1Boolean toleranceControlled, fmi1Real relativeTolerance, fmi1EventInfo* eventInfo);
typedef fmi1Status    fmi1GetDerivativesTYPE             (fmi1Component c, fmi1Real derivatives[], size_t nx);
typedef fmi1Status    fmi1GetEventIndicatorsTYPE         (fmi1Component c, fmi1Real eventIndicators[], size_t ni);
typedef fmi1Status    fmi1EventUpdateTYPE                (fmi1Component c, fmi1Boolean intermediateResults, fmi1EventInfo* eventInfo);
typedef fmi1Status    fmi1GetContinuousStatesTYPE        (fmi1Component c, fmi1Real states[], size_t nx);
typedef fmi1Status    fmi1GetNominalContinuousStatesTYPE (fmi1Component c, fmi1Real x_nominal[], size_t nx);
typedef fmi1Status    fmi1GetStateValueReferencesTYPE    (fmi1Component c, fmi1ValueReference vrx[], size_t nx);
typedef fmi1Status    fmi1TerminateTYPE                  (fmi1Component c);

/***************************************************
 FMI 1.0 for Co-Simulation Functions
****************************************************/
typedef const char*   fmi1GetTypesPlatformTYPE         ();
typedef fmi1Component fmi1InstantiateSlaveTYPE         (fmi1String  instanceName, fmi1String  fmuGUID, fmi1String  fmuLocation,fmi1String  mimeType, fmi1Real timeout, fmi1Boolean visible, fmi1Boolean interactive, fmi1CallbackFunctions functions, fmi1Boolean loggingOn);
typedef fmi1Status    fmi1InitializeSlaveTYPE          (fmi1Component c, fmi1Real tStart, fmi1Boolean StopTimeDefined, fmi1Real tStop);
typedef fmi1Status    fmi1TerminateSlaveTYPE           (fmi1Component c);
typedef fmi1Status    fmi1ResetSlaveTYPE               (fmi1Component c);
typedef void          fmi1FreeSlaveInstanceTYPE        (fmi1Component c);
typedef fmi1Status    fmi1SetRealInputDerivativesTYPE  (fmi1Component c, const  fmi1ValueReference vr[], size_t nvr, const fmi1Integer order[], const  fmi1Real value[]);
typedef fmi1Status    fmi1GetRealOutputDerivativesTYPE (fmi1Component c, const fmi1ValueReference vr[], size_t  nvr, const fmi1Integer order[], fmi1Real value[]);
typedef fmi1Status    fmi1CancelStepTYPE               (fmi1Component c);
typedef fmi1Status    fmi1DoStepTYPE                   (fmi1Component c, fmi1Real currentCommunicationPoint, fmi1Real communicationStepSize, fmi1Boolean newStep);
typedef fmi1Status    fmi1GetStatusTYPE                (fmi1Component c, const fmi1StatusKind s, fmi1Status*  value);
typedef fmi1Status    fmi1GetRealStatusTYPE            (fmi1Component c, const fmi1StatusKind s, fmi1Real*    value);
typedef fmi1Status    fmi1GetIntegerStatusTYPE         (fmi1Component c, const fmi1StatusKind s, fmi1Integer* value);
typedef fmi1Status    fmi1GetBooleanStatusTYPE         (fmi1Component c, const fmi1StatusKind s, fmi1Boolean* value);
typedef fmi1Status    fmi1GetStringStatusTYPE          (fmi1Component c, const fmi1StatusKind s, fmi1String*  value);

#endif // FMI1FUNCTIONS_H
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <Windows.h>
#endif

#include <stdbool.h>
#include <stddef.h>

#ifndef FMI_MAX_MESSAGE_LENGTH
#define FMI_MAX_MESSAGE_LENGTH 4096
#endif

#ifndef FMI_STATIC
#define FMI_STATIC
#endif

typedef enum {
    FMIOK,
    FMIWarning,
    FMIDiscard,
    FMIError,
    FMIFatal,
    FMIPending
} FMIStatus;

typedef enum {

    // FMI 3.0 variable types
    FMIFloat32Type,
    FMIDiscreteFloat32Type,
    FMIFloat64Type,
    FMIDiscreteFloat64Type,
    FMIInt8Type,
    FMIUInt8Type,
    FMIInt16Type,
    FMIUInt16Type,
    FMIInt32Type,
    FMIUInt32Type,
    FMIInt64Type,
    FMIUInt64Type,
    FMIBooleanType,
    FMIStringType,
    FMIBinaryType,
    FMIClockType,

    // non-variable types used as parameters
    FMIValueReferenceType,
    FMISizeTType,

    // aliases for FMI 1.0 and 2.0 variable types
    FMIRealType         = FMIFloat64Type,
    FMIDiscreteRealType = FMIDiscreteFloat64Type,
    FMIIntegerType      = FMIInt32Type,

} FMIVariableType;

typedef enum {
    FMIVersion1 = 1,
    FMIVersion2 = 2,
    FMIVersion3 = 3
} FMIVersion;

typedef enum {
    FMIModelExchange,
    FMICoSimulation,
    FMIScheduledExecution
} FMIInterfaceType;

typedef enum {
    FMI2StartAndEndState        = 1 << 0,
    FMI2InstantiatedState       = 1 << 1,
    FMI2InitializationModeState = 1 << 2,

    // model exchange states
    FMI2EventModeState          = 1 << 3,
    FMI2ContinuousTimeModeState = 1 << 4,

    // co-simulation states
    FMI2StepCompleteState       = 1 << 5,
    FMI2StepInProgressState     = 1 << 6,
    FMI2StepFailedState         = 1 << 7,
    FMI2StepCanceledState       = 1 << 8,

    FMI2TerminatedState         = 1 << 9,
    FMI2ErrorState              = 1 << 10,
    FMI2FatalState              = 1 << 11,
} FMI2State;

typedef unsigned int FMIValueReference;

typedef struct FMIInstance_ FMIInstance;

typedef struct FMI1Functions_ FMI1Functions;

typedef struct FMI2Functions_ FMI2Functions;

typedef struct FMI3Functions_ FMI3Functions;

typedef void FMILogFunctionCall(FMIInstance *instance, FMIStatus status, const char *message);

typedef void FMILogMessage(FMIInstance *instance, FMIStatus status, const char *category, const char *message);

struct FMIInstance_ {

    FMI1Functions *fmi1Functions;
    FMI2Functions *fmi2Functions;
    FMI3Functions *fmi3Functions;

#ifdef _WIN32
    HMODULE libraryHandle;
#else
    void *libraryHandle;
#endif

    void *userData;

    FMILogMessage      *logMessage;
    FMILogFunctionCall *logFunctionCall;

    double time;

    char* logMessageBuffer;
    size_t logMessageBufferSize;
    size_t logMessageBufferPosition;

    void *component;

    const char *name;

    bool logFMICalls;

    FMI2State state;

    FMIStatus status;

    FMIVersion fmiVersion;

    FMIInterfaceType interfaceType;

};

FMI_STATIC FMIInstance* FMICreateInstance(const char* instanceName, const char* libraryPath, FMILogMessage* logMessage, FMILogFunctionCall* logFunctionCall);

FMI_STATIC void FMIFreeInstance(FMIInstance *instance);

FMI_STATIC void FMIClearLogMessageBuffer(FMIInstance* instance);

FMI_STATIC void FMIAppendToLogMessageBuffer(FMIInstance* instance, const char* format, ...);

FMI_STATIC void FMIAppendArrayToLogMessageBuffer(FMIInstance* instance, const void* values, size_t nValues, const size_t sizes[], FMIVariableType variableType);

FMI_STATIC FMIStatus FMIURIToPath(const char *uri, char *path, const size_t pathLength);

FMI_STATIC FMIStatus FMIPathToURI(const char *path, char *uri, const size_t uriLength);

FMI_STATIC FMIStatus FMIPlatformBinaryPath(const char *unzipdir, const char *modelIdentifier, FMIVersion fmiVersion, char *platformBinaryPath, size_t size);

#ifdef __cplusplus
}  /* end of extern "C" { */
#endif

struct FMI1Functions_ {

    /***************************************************
     Common Functions for FMI 1.0
    ****************************************************/

    fmi1CallbackFunctions   callbacks;
    fmi1EventInfo           eventInfo;

    fmi1SetRealTYPE         *fmi1SetReal;
    fmi1SetIntegerTYPE      *fmi1SetInteger;
    fmi1SetBooleanTYPE      *fmi1SetBoolean;
    fmi1SetStringTYPE       *fmi1SetString;
    fmi1GetRealTYPE         *fmi1GetReal;
    fmi1GetIntegerTYPE      *fmi1GetInteger;
    fmi1GetBooleanTYPE      *fmi1GetBoolean;
    fmi1GetStringTYPE       *fmi1GetString;
    fmi1SetDebugLoggingTYPE *fmi1SetDebugLogging;

    /***************************************************
     FMI 1.0 for Model Exchange Functions
    ****************************************************/

    fmi1GetModelTypesPlatformTYPE      *fmi1GetModelTypesPlatform;
    fmi1GetVersionTYPE                 *fmi1GetVersion;
    fmi1InstantiateModelTYPE           *fmi1InstantiateModel;
    fmi1FreeModelInstanceTYPE          *fmi1FreeModelInstance;
    fmi1SetTimeTYPE                    *fmi1SetTime;
    fmi1SetContinuousStatesTYPE        *fmi1SetContinuousStates;
    fmi1CompletedIntegratorStepTYPE    *fmi1CompletedIntegratorStep;
    fmi1InitializeTYPE                 *fmi1Initialize;
    fmi1GetDerivativesTYPE             *fmi1GetDerivatives;
    fmi1GetEventIndicatorsTYPE         *fmi1GetEventIndicators;
    fmi1EventUpdateTYPE                *fmi1EventUpdate;
    fmi1GetContinuousStatesTYPE        *fmi1GetContinuousStates;
    fmi1GetNominalContinuousStatesTYPE *fmi1GetNominalContinuousStates;
    fmi1GetStateValueReferencesTYPE    *fmi1GetStateValueReferences;
    fmi1TerminateTYPE                  *fmi1Terminate;

    /***************************************************
     FMI 1.0 for Co-Simulation Functions
    ****************************************************/

    fmi1GetTypesPlatformTYPE         *fmi1GetTypesPlatform;
    fmi1InstantiateSlaveTYPE         *fmi1InstantiateSlave;
    fmi1InitializeSlaveTYPE          *fmi1InitializeSlave;
    fmi1TerminateSlaveTYPE           *fmi1TerminateSlave;
    fmi1ResetSlaveTYPE               *fmi1ResetSlave;
    fmi1FreeSlaveInstanceTYPE        *fmi1FreeSlaveInstance;
    fmi1SetRealInputDerivativesTYPE  *fmi1SetRealInputDerivatives;
    fmi1GetRealOutputDerivativesTYPE *fmi1GetRealOutputDerivatives;
    fmi1CancelStepTYPE               *fmi1CancelStep;
    fmi1DoStepTYPE                   *fmi1DoStep;
    fmi1GetStatusTYPE                *fmi1GetStatus;
    fmi1GetRealStatusTYPE            *fmi1GetRealStatus;
    fmi1GetIntegerStatusTYPE         *fmi1GetIntegerStatus;
    fmi1GetBooleanStatusTYPE         *fmi1GetBooleanStatus;
    fmi1GetStringStatusTYPE          *fmi1GetStringStatus;

};

/***************************************************
 Common Functions for FMI 1.0
****************************************************/
FMI_STATIC FMIStatus    FMI1SetReal         (FMIInstance *instance, const fmi1ValueReference vr[], size_t nvr, const fmi1Real    value[]);
FMI_STATIC FMIStatus    FMI1SetInteger      (FMIInstance *instance, const fmi1ValueReference vr[], size_t nvr, const fmi1Integer value[]);
FMI_STATIC FMIStatus    FMI1SetBoolean      (FMIInstance *instance, const fmi1ValueReference vr[], size_t nvr, const fmi1Boolean value[]);
FMI_STATIC FMIStatus    FMI1SetString       (FMIInstance *instance, const fmi1ValueReference vr[], size_t nvr, const fmi1String  value[]);
FMI_STATIC FMIStatus    FMI1GetReal         (FMIInstance *instance, const fmi1ValueReference vr[], size_t nvr,       fmi1Real    value[]);
FMI_STATIC FMIStatus    FMI1GetInteger      (FMIInstance *instance, const fmi1ValueReference vr[], size_t nvr,       fmi1Integer value[]);
FMI_STATIC FMIStatus    FMI1GetBoolean      (FMIInstance *instance, const fmi1ValueReference vr[], size_t nvr,       fmi1Boolean value[]);
FMI_STATIC FMIStatus    FMI1GetString       (FMIInstance *instance, const fmi1ValueReference vr[], size_t nvr,       fmi1String  value[]);
FMI_STATIC FMIStatus    FMI1SetDebugLogging (FMIInstance *instance, fmi1Boolean loggingOn);

/***************************************************
 FMI 1.0 for Model Exchange Functions
****************************************************/
FMI_STATIC const char*  FMI1GetModelTypesPlatform      (FMIInstance *instance);
FMI_STATIC const char*  FMI1GetVersion                 (FMIInstance *instance);
FMI_STATIC FMIStatus    FMI1InstantiateModel           (FMIInstance *instance, fmi1String modelIdentifier, fmi1String GUID, fmi1Boolean loggingOn);
FMI_STATIC void         FMI1FreeModelInstance          (FMIInstance *instance);
FMI_STATIC FMIStatus    FMI1SetTime                    (FMIInstance *instance, fmi1Real time);
FMI_STATIC FMIStatus    FMI1SetContinuousStates        (FMIInstance *instance, const fmi1Real x[], size_t nx);
FMI_STATIC FMIStatus    FMI1CompletedIntegratorStep    (FMIInstance *instance, fmi1Boolean* callEventUpdate);
FMI_STATIC FMIStatus    FMI1Initialize                 (FMIInstance *instance, fmi1Boolean toleranceControlled, fmi1Real relativeTolerance, fmi1EventInfo* eventInfo);
FMI_STATIC FMIStatus    FMI1GetDerivatives             (FMIInstance *instance, fmi1Real derivatives[], size_t nx);
FMI_STATIC FMIStatus    FMI1GetEventIndicators         (FMIInstance *instance, fmi1Real eventIndicators[], size_t ni);
FMI_STATIC FMIStatus    FMI1EventUpdate                (FMIInstance *instance, fmi1Boolean intermediateResults, fmi1EventInfo* eventInfo);
FMI_STATIC FMIStatus    FMI1GetContinuousStates        (FMIInstance *instance, fmi1Real states[], size_t nx);
FMI_STATIC FMIStatus    FMI1GetNominalContinuousStates (FMIInstance *instance, fmi1Real x_nominal[], size_t nx);
FMI_STATIC FMIStatus    FMI1GetStateValueReferences    (FMIInstance *instance, fmi1ValueReference vrx[], size_t nx);
FMI_STATIC FMIStatus    FMI1Terminate                  (FMIInstance *instance);

/***************************************************
 FMI 1.0 for Co-Simulation Functions
****************************************************/
FMI_STATIC const char*  FMI1GetTypesPlatform         (FMIInstance *instance);
FMI_STATIC FMIStatus    FMI1InstantiateSlave         (FMIInstance *instance, fmi1String modelIdentifier, fmi1String fmuGUID, fmi1String fmuLocation, fmi1String  mimeType, fmi1Real timeout, fmi1Boolean visible, fmi1Boolean interactive, fmi1Boolean loggingOn);
FMI_STATIC FMIStatus    FMI1InitializeSlave          (FMIInstance *instance, fmi1Real tStart, fmi1Boolean StopTimeDefined, fmi1Real tStop);
FMI_STATIC FMIStatus    FMI1TerminateSlave           (FMIInstance *instance);
FMI_STATIC FMIStatus    FMI1ResetSlave               (FMIInstance *instance);
FMI_STATIC void         FMI1FreeSlaveInstance        (FMIInstance *instance);
FMI_STATIC FMIStatus    FMI1SetRealInputDerivatives  (FMIInstance *instance, const fmi1ValueReference vr[], size_t nvr, const fmi1Integer order[], const fmi1Real value[]);
FMI_STATIC FMIStatus    FMI1GetRealOutputDerivatives (FMIInstance *instance, const fmi1ValueReference vr[], size_t nvr, const fmi1Integer order[],       fmi1Real value[]);
FMI_STATIC FMIStatus    FMI1CancelStep               (FMIInstance *instance);
FMI_STATIC FMIStatus    FMI1DoStep                   (FMIInstance *instance, fmi1Real currentCommunicationPoint, fmi1Real communicationStepSize, fmi1Boolean newStep);
FMI_STATIC FMIStatus    FMI1GetStatus                (FMIInstance *instance, const fmi1StatusKind s, fmi1Status*  value);
FMI_STATIC FMIStatus    FMI1GetRealStatus            (FMIInstance *instance, const fmi1StatusKind s, fmi1Real*    value);
FMI_STATIC FMIStatus    FMI1GetIntegerStatus         (FMIInstance *instance, const fmi1StatusKind s, fmi1Integer* value);
FMI_STATIC FMIStatus    FMI1GetBooleanStatus         (FMIInstance *instance, const fmi1StatusKind s, fmi1Boolean* value);
FMI_STATIC FMIStatus    FMI1GetStringStatus          (FMIInstance *instance, const fmi1StatusKind s, fmi1String*  value);

#ifdef __cplusplus
}  /* end of extern "C" { */
#endif


// keep track of the current instance for the log callback
static FMIInstance *currentInstance = NULL;

static void cb_logMessage1(fmi1Component c, fmi1String instanceName, fmi1Status status, fmi1String category, fmi1String message, ...) {

    if (!currentInstance) return;

    char buf[FMI_MAX_MESSAGE_LENGTH];

    va_list args;

    va_start(args, message);
    vsnprintf(buf, FMI_MAX_MESSAGE_LENGTH, message, args);
    va_end(args);

    if (!currentInstance->logMessage) return;

    currentInstance->logMessage(currentInstance, (FMIStatus)status, category, buf);
}

#define MAX_SYMBOL_LENGTH 256

static void *loadSymbol(FMIInstance *instance, const char *prefix, const char *name) {
    char fname[MAX_SYMBOL_LENGTH];
    strcpy(fname, prefix);
    strcat(fname, name);
#ifdef _WIN32
    void *addr = GetProcAddress(instance->libraryHandle, fname);
#else
    void *addr = dlsym(instance->libraryHandle, fname);
#endif
    if (!addr) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "Failed to load function \"%s\".", fname);
        instance->logFunctionCall(instance, FMIError, instance->logMessageBuffer);
    }
    return addr;
}

#define LOAD_SYMBOL(f) \
do { \
    instance->fmi1Functions->fmi1 ## f = (fmi1 ## f ## TYPE*)loadSymbol(instance, modelIdentifier, "_fmi" #f); \
    if (!instance->fmi1Functions->fmi1 ## f) { \
        status = FMIError; \
        goto fail; \
    } \
} while (0)

#define CALL(f) \
do { \
    currentInstance = instance; \
    FMIStatus status = (FMIStatus)instance->fmi1Functions->fmi1 ## f (instance->component); \
    if (instance->logFunctionCall) { \
        instance->logFunctionCall(instance, status, "fmi" #f "()"); \
    } \
    return status; \
} while (0)

#define CALL_ARGS(f, m, ...) \
do { \
    currentInstance = instance; \
    FMIStatus status = (FMIStatus)instance->fmi1Functions->fmi1 ## f (instance->component, __VA_ARGS__); \
    if (instance->logFunctionCall) { \
        FMIClearLogMessageBuffer(instance); \
        FMIAppendToLogMessageBuffer(instance, "fmi" #f "(" m ")", __VA_ARGS__); \
        instance->logFunctionCall(instance, status, instance->logMessageBuffer); \
    } \
    return status; \
} while (0)

#define CALL_ARRAY(s, t) \
do { \
    currentInstance = instance; \
    FMIStatus status = (FMIStatus)instance->fmi1Functions->fmi1 ## s ## t(instance->component, vr, nvr, value); \
    if (instance->logFunctionCall) { \
        FMIClearLogMessageBuffer(instance); \
        FMIAppendToLogMessageBuffer(instance, "fmi" #s #t "(vr={"); \
        FMIAppendArrayToLogMessageBuffer(instance, vr, nvr, NULL, FMIValueReferenceType); \
        FMIAppendToLogMessageBuffer(instance, "}, nvr=%zu, value={", nvr); \
        FMIAppendArrayToLogMessageBuffer(instance, value, nvr, NULL, FMI ## t ## Type); \
        FMIAppendToLogMessageBuffer(instance, "})"); \
        instance->logFunctionCall(instance, status, instance->logMessageBuffer); \
    } \
    return status; \
} while (0)

/***************************************************
 Common Functions for FMI 1.0
****************************************************/

FMIStatus FMI1SetReal(FMIInstance *instance, const fmi1ValueReference vr[], size_t nvr, const fmi1Real    value[]) {
    CALL_ARRAY(Set, Real);
}

FMIStatus FMI1SetInteger(FMIInstance *instance, const fmi1ValueReference vr[], size_t nvr, const fmi1Integer value[]) {
    CALL_ARRAY(Set, Integer);
}

FMIStatus FMI1SetBoolean(FMIInstance *instance, const fmi1ValueReference vr[], size_t nvr, const fmi1Boolean value[]) {
    CALL_ARRAY(Set, Boolean);
}

FMIStatus FMI1SetString(FMIInstance *instance, const fmi1ValueReference vr[], size_t nvr, const fmi1String  value[]) {
    CALL_ARRAY(Set, String);
}

FMIStatus FMI1GetReal(FMIInstance *instance, const fmi1ValueReference vr[], size_t nvr, fmi1Real    value[]) {
    CALL_ARRAY(Get, Real);
}

FMIStatus FMI1GetInteger(FMIInstance *instance, const fmi1ValueReference vr[], size_t nvr, fmi1Integer value[]) {
    CALL_ARRAY(Get, Integer);
}

FMIStatus FMI1GetBoolean(FMIInstance *instance, const fmi1ValueReference vr[], size_t nvr, fmi1Boolean value[]) {
    CALL_ARRAY(Get, Boolean);
}

FMIStatus FMI1GetString(FMIInstance *instance, const fmi1ValueReference vr[], size_t nvr, fmi1String  value[]) {
    CALL_ARRAY(Get, String);
}

FMIStatus FMI1SetDebugLogging(FMIInstance *instance, fmi1Boolean loggingOn) {
    CALL_ARGS(SetDebugLogging, "loggingOn=%d", loggingOn);
}


/***************************************************
 FMI 1.0 for Model Exchange Functions
****************************************************/

const char* FMI1GetModelTypesPlatform(FMIInstance *instance) {
    currentInstance = instance;
    if (instance->logFunctionCall) {
        instance->logFunctionCall(instance, FMIOK, "fmiGetModelTypesPlatform()");
    }
    return instance->fmi1Functions->fmi1GetModelTypesPlatform();
}

const char* FMI1GetVersion(FMIInstance *instance) {
    currentInstance = instance;
    if (instance->logFunctionCall) {
        instance->logFunctionCall(instance, FMIOK, "fmiGetVersion()");
    }
    return instance->fmi1Functions->fmi1GetVersion();
}

FMIStatus FMI1InstantiateModel(FMIInstance *instance, fmi1String modelIdentifier, fmi1String GUID, fmi1Boolean loggingOn) {

    currentInstance = instance;

    instance->fmiVersion = FMIVersion1;

    instance->interfaceType = FMIModelExchange;

    FMIStatus status = FMIOK;

    instance->fmi1Functions = calloc(1, sizeof(FMI1Functions));

    if (!instance->fmi1Functions) {
        return FMIError;
    }

    /***************************************************
     Common Functions for FMI 1.0
    ****************************************************/
    LOAD_SYMBOL(SetReal);
    LOAD_SYMBOL(SetInteger);
    LOAD_SYMBOL(SetBoolean);
    LOAD_SYMBOL(SetString);
    LOAD_SYMBOL(GetReal);
    LOAD_SYMBOL(GetInteger);
    LOAD_SYMBOL(GetBoolean);
    LOAD_SYMBOL(GetString);
    LOAD_SYMBOL(SetDebugLogging);

    /***************************************************
        FMI 1.0 for Model Exchange Functions
    ****************************************************/
    LOAD_SYMBOL(GetModelTypesPlatform);
    LOAD_SYMBOL(GetVersion);
    LOAD_SYMBOL(InstantiateModel);
    LOAD_SYMBOL(FreeModelInstance);
    LOAD_SYMBOL(SetTime);
    LOAD_SYMBOL(SetContinuousStates);
    LOAD_SYMBOL(CompletedIntegratorStep);
    LOAD_SYMBOL(Initialize);
    LOAD_SYMBOL(GetDerivatives);
    LOAD_SYMBOL(GetEventIndicators);
    LOAD_SYMBOL(EventUpdate);
    LOAD_SYMBOL(GetContinuousStates);
    LOAD_SYMBOL(GetNominalContinuousStates);
    LOAD_SYMBOL(GetStateValueReferences);
    LOAD_SYMBOL(Terminate);

    instance->fmi1Functions->callbacks.logger         = cb_logMessage1;
    instance->fmi1Functions->callbacks.allocateMemory = calloc;
    instance->fmi1Functions->callbacks.freeMemory     = free;
    instance->fmi1Functions->callbacks.stepFinished   = NULL;

    instance->component = instance->fmi1Functions->fmi1InstantiateModel(instance->name, GUID, instance->fmi1Functions->callbacks, loggingOn);

    status = instance->component ? FMIOK : FMIError;

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance,
            "fmiInstantiateModel(instanceName=\"%s\", GUID=\"%s\", functions=0x%p, loggingOn=%d)",
            instance->name, GUID, &instance->fmi1Functions->callbacks, loggingOn);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

fail:
    return status;
}

void FMI1FreeModelInstance(FMIInstance *instance) {

    if (!instance) {
        return;
    }

    currentInstance = instance;

    instance->fmi1Functions->fmi1FreeModelInstance(instance->component);

    if (instance->logFunctionCall) {
        instance->logFunctionCall(instance, FMIOK, "fmiFreeModelInstance()");
    }
}

FMIStatus FMI1SetTime(FMIInstance *instance, fmi1Real time) {
    CALL_ARGS(SetTime, "time=%.16g", time);
}

FMIStatus FMI1SetContinuousStates(FMIInstance *instance, const fmi1Real x[], size_t nx) {

    currentInstance = instance;

    const FMIStatus status = (FMIStatus)instance->fmi1Functions->fmi1SetContinuousStates(instance->component, x, nx);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmiSetContinuousStates(x={");
        FMIAppendArrayToLogMessageBuffer(instance, x, nx, NULL, FMIRealType);
        FMIAppendToLogMessageBuffer(instance, "}, nx=%zu)", nx);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI1CompletedIntegratorStep(FMIInstance *instance, fmi1Boolean* callEventUpdate) {

    currentInstance = instance;

    const FMIStatus status = (FMIStatus)instance->fmi1Functions->fmi1CompletedIntegratorStep(instance->component, callEventUpdate);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmiCompletedIntegratorStep(callEventUpdate=%d)", *callEventUpdate);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI1Initialize(FMIInstance *instance, fmi1Boolean toleranceControlled, fmi1Real relativeTolerance, fmi1EventInfo* eventInfo) {

    currentInstance = instance;

    const FMIStatus status = (FMIStatus)instance->fmi1Functions->fmi1Initialize(instance->component, toleranceControlled, relativeTolerance, eventInfo);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance,
            "fmiInitialize(toleranceControlled=%d, relativeTolerance=%.16g, eventInfo={iterationConverged=%d, stateValueReferencesChanged=%d, stateValuesChanged=%d, terminateSimulation=%d, upcomingTimeEvent=%d, nextEventTime=%.16g})",
            toleranceControlled, relativeTolerance, eventInfo->iterationConverged, eventInfo->stateValueReferencesChanged, eventInfo->stateValuesChanged, eventInfo->terminateSimulation, eventInfo->upcomingTimeEvent, eventInfo->nextEventTime);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI1GetDerivatives(FMIInstance *instance, fmi1Real derivatives[], size_t nx) {

    currentInstance = instance;

    const FMIStatus status = (FMIStatus)instance->fmi1Functions->fmi1GetDerivatives(instance->component, derivatives, nx);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmiGetDerivatives(derivatives={");
        FMIAppendArrayToLogMessageBuffer(instance, derivatives, nx, NULL, FMIRealType);
        FMIAppendToLogMessageBuffer(instance, "}, nx=%zu)", nx);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI1GetEventIndicators(FMIInstance *instance, fmi1Real eventIndicators[], size_t ni) {

    currentInstance = instance;

    const FMIStatus status = (FMIStatus)instance->fmi1Functions->fmi1GetEventIndicators(instance->component, eventIndicators, ni);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmiGetEventIndicators(eventIndicators={");
        FMIAppendArrayToLogMessageBuffer(instance, eventIndicators, ni, NULL, FMIRealType);
        FMIAppendToLogMessageBuffer(instance, "}, ni=%zu)", ni);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI1EventUpdate(FMIInstance *instance, fmi1Boolean intermediateResults, fmi1EventInfo* eventInfo) {

    currentInstance = instance;

    const FMIStatus status = (FMIStatus)instance->fmi1Functions->fmi1EventUpdate(instance->component, intermediateResults, eventInfo);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance,
            "fmiEventUpdate(intermediateResults=%d, eventInfo={iterationConverged=%d, stateValueReferencesChanged=%d, stateValuesChanged=%d, terminateSimulation=%d, upcomingTimeEvent=%d, nextEventTime=%.16g})",
            intermediateResults, eventInfo->iterationConverged, eventInfo->stateValueReferencesChanged, eventInfo->stateValuesChanged, eventInfo->terminateSimulation, eventInfo->upcomingTimeEvent, eventInfo->nextEventTime);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI1GetContinuousStates(FMIInstance *instance, fmi1Real states[], size_t nx) {

    currentInstance = instance;

    const FMIStatus status = (FMIStatus)instance->fmi1Functions->fmi1GetContinuousStates(instance->component, states, nx);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmiGetContinuousStates(states={");
        FMIAppendArrayToLogMessageBuffer(instance, states, nx, NULL, FMIRealType);
        FMIAppendToLogMessageBuffer(instance, "}, nx=%zu)", nx);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI1GetNominalContinuousStates(FMIInstance *instance, fmi1Real x_nominal[], size_t nx) {

    currentInstance = instance;

    const FMIStatus status = (FMIStatus)instance->fmi1Functions->fmi1GetNominalContinuousStates(instance->component, x_nominal, nx);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmiGetNominalContinuousStates(x_nominal={");
        FMIAppendArrayToLogMessageBuffer(instance, x_nominal, nx, NULL, FMIRealType);
        FMIAppendToLogMessageBuffer(instance, "}, nx=%zu)", nx);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI1GetStateValueReferences(FMIInstance *instance, fmi1ValueReference vrx[], size_t nx) {

    currentInstance = instance;

    const FMIStatus status = (FMIStatus)instance->fmi1Functions->fmi1GetStateValueReferences(instance->component, vrx, nx);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmiGetStateValueReferences(vrx={");
        FMIAppendArrayToLogMessageBuffer(instance, vrx, nx, NULL, FMIValueReferenceType);
        FMIAppendToLogMessageBuffer(instance, "}, nx=%zu)", nx);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI1Terminate(FMIInstance *instance) {
    CALL(Terminate);
}

/***************************************************
 FMI 1.0 for Co-Simulation Functions
****************************************************/

const char*   FMI1GetTypesPlatform(FMIInstance *instance) {

    currentInstance = instance;

    if (instance->logFunctionCall) {
        instance->logFunctionCall(instance, FMIOK, "fmiGetTypesPlatform()");
    }

    return instance->fmi1Functions->fmi1GetTypesPlatform();
}

FMIStatus FMI1InstantiateSlave(FMIInstance *instance, fmi1String modelIdentifier, fmi1String fmuGUID, fmi1String fmuLocation, fmi1String  mimeType, fmi1Real timeout, fmi1Boolean visible, fmi1Boolean interactive, fmi1Boolean loggingOn) {

    currentInstance = instance;

    instance->fmiVersion = FMIVersion1;

    instance->interfaceType = FMICoSimulation;

    FMIStatus status = FMIOK;

    instance->fmi1Functions = calloc(1, sizeof(FMI1Functions));

    if (!instance->fmi1Functions) {
        return FMIError;
    }

    /***************************************************
     Common Functions for FMI 1.0
    ****************************************************/
    LOAD_SYMBOL(SetReal);
    LOAD_SYMBOL(SetInteger);
    LOAD_SYMBOL(SetBoolean);
    LOAD_SYMBOL(SetString);
    LOAD_SYMBOL(GetReal);
    LOAD_SYMBOL(GetInteger);
    LOAD_SYMBOL(GetBoolean);
    LOAD_SYMBOL(GetString);
    LOAD_SYMBOL(SetDebugLogging);

    /***************************************************
     FMI 1.0 for Co-Simulation Functions
    ****************************************************/
    LOAD_SYMBOL(GetTypesPlatform);
    LOAD_SYMBOL(InstantiateSlave);
    LOAD_SYMBOL(InitializeSlave);
    LOAD_SYMBOL(TerminateSlave);
    LOAD_SYMBOL(ResetSlave);
    LOAD_SYMBOL(FreeSlaveInstance);
    LOAD_SYMBOL(SetRealInputDerivatives);
    LOAD_SYMBOL(GetRealOutputDerivatives);
    LOAD_SYMBOL(CancelStep);
    LOAD_SYMBOL(DoStep);
    LOAD_SYMBOL(GetStatus);
    LOAD_SYMBOL(GetRealStatus);
    LOAD_SYMBOL(GetIntegerStatus);
    LOAD_SYMBOL(GetBooleanStatus);
    LOAD_SYMBOL(GetStringStatus);

    instance->fmi1Functions->callbacks.logger         = cb_logMessage1;
    instance->fmi1Functions->callbacks.allocateMemory = calloc;
    instance->fmi1Functions->callbacks.freeMemory     = free;
    instance->fmi1Functions->callbacks.stepFinished   = NULL;

    instance->component = instance->fmi1Functions->fmi1InstantiateSlave(instance->name, fmuGUID, fmuLocation, mimeType, timeout, visible, interactive, instance->fmi1Functions->callbacks, loggingOn);

    status = instance->component ? FMIOK : FMIError;

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance,
            "fmiInstantiateSlave(instanceName=\"%s\", fmuGUID=\"%s\", fmuLocation=\"%s\", mimeType=\"%s\", timeout=%.16g, visible=%d, interactive=%d, functions=0x%p, loggingOn=%d)",
            instance->name, fmuGUID, fmuLocation, mimeType, timeout, visible, interactive, &instance->fmi1Functions->callbacks, loggingOn);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

fail:
    return status;
}

FMIStatus FMI1InitializeSlave(FMIInstance *instance, fmi1Real tStart, fmi1Boolean stopTimeDefined, fmi1Real tStop) {
    CALL_ARGS(InitializeSlave, "tStart=%.16g, stopTimeDefined=%d, tStop=%.16g", tStart, stopTimeDefined, tStop);
}

FMIStatus FMI1TerminateSlave(FMIInstance *instance) {
    CALL(TerminateSlave);
}

FMIStatus FMI1ResetSlave(FMIInstance *instance) {
    CALL(ResetSlave);
}

void FMI1FreeSlaveInstance(FMIInstance *instance) {

    if (!instance) {
        return;
    }

    currentInstance = instance;

    instance->fmi1Functions->fmi1FreeSlaveInstance(instance->component);

    if (instance->logFunctionCall) {
        instance->logFunctionCall(instance, FMIOK, "fmiFreeSlaveInstance()");
    }
}

FMIStatus FMI1SetRealInputDerivatives(FMIInstance *instance, const fmi1ValueReference vr[], size_t nvr, const fmi1Integer order[], const fmi1Real value[]) {
    CALL_ARGS(SetRealInputDerivatives, "vr=0x%p, nvr=%zu, order=0x%p, value=0x%p", vr, nvr, order, value);
}

FMIStatus FMI1GetRealOutputDerivatives(FMIInstance *instance, const fmi1ValueReference vr[], size_t nvr, const fmi1Integer order[], fmi1Real value[]) {
    CALL_ARGS(GetRealOutputDerivatives, "vr=0x%p, nvr=%zu, order=0x%p, value=0x%p", vr, nvr, order, value);
}

FMIStatus FMI1CancelStep(FMIInstance *instance) {
    CALL(CancelStep);
}

FMIStatus FMI1DoStep(FMIInstance *instance, fmi1Real currentCommunicationPoint, fmi1Real communicationStepSize, fmi1Boolean newStep) {

    currentInstance = instance;

    instance->time = currentCommunicationPoint + communicationStepSize;

    CALL_ARGS(DoStep, "currentCommunicationPoint=%.16g, communicationStepSize=%.16g, newStep=%d",
        currentCommunicationPoint, communicationStepSize, newStep);
}

FMIStatus FMI1GetStatus(FMIInstance *instance, const fmi1StatusKind s, fmi1Status* value) {
    CALL_ARGS(GetStatus, "s=%d, value=0x%p", s, value);
}

FMIStatus FMI1GetRealStatus(FMIInstance *instance, const fmi1StatusKind s, fmi1Real* value) {
    CALL_ARGS(GetRealStatus, "s=%d, value=0x%p", s, value);
}

FMIStatus FMI1GetIntegerStatus(FMIInstance *instance, const fmi1StatusKind s, fmi1Integer* value) {
    CALL_ARGS(GetIntegerStatus, "s=%d, value=0x%p", s, value);
}

FMIStatus FMI1GetBooleanStatus(FMIInstance *instance, const fmi1StatusKind s, fmi1Boolean* value) {
    CALL_ARGS(GetBooleanStatus, "s=%d, value=0x%p", s, value);
}

FMIStatus FMI1GetStringStatus(FMIInstance *instance, const fmi1StatusKind s, fmi1String* value) {
    CALL_ARGS(GetStringStatus, "s=%d, value=0x%p", s, value);
}

#undef LOAD_SYMBOL
#undef CALL
#undef CALL_ARGS
#undef CALL_ARRAY
#ifdef _WIN32
#include <direct.h>
#include "Shlwapi.h"
#pragma comment(lib, "shlwapi.lib")
#define strdup _strdup
#define INTERNET_MAX_URL_LENGTH 2083 // from wininet.h
#else
#include <stdarg.h>
#include <dlfcn.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#ifndef fmi2Functions_h
#define fmi2Functions_h

/* This header file must be utilized when compiling a FMU.
   It defines all functions of the
         FMI 2.0.1 Model Exchange and Co-Simulation Interface.

   In order to have unique function names even if several FMUs
   are compiled together (e.g. for embedded systems), every "real" function name
   is constructed by prepending the function name by "FMI2_FUNCTION_PREFIX".
   Therefore, the typical usage is:

      #define FMI2_FUNCTION_PREFIX MyModel_
      #include "fmi2Functions.h"

   As a result, a function that is defined as "fmi2GetDerivatives" in this header file,
   is actually getting the name "MyModel_fmi2GetDerivatives".

   This only holds if the FMU is shipped in C source code, or is compiled in a
   static link library. For FMUs compiled in a DLL/sharedObject, the "actual" function
   names are used and "FMI2_FUNCTION_PREFIX" must not be defined.

   Revisions:
   - Sep. 29, 2019: License changed to 2-clause BSD License (without extensions)
   - Apr.  9, 2014: All prefixes "fmi" renamed to "fmi2" (decision from April 8)
   - Mar. 26, 2014: FMI_Export set to empty value if FMI_Export and FMI_FUNCTION_PREFIX
                    are not defined (#173)
   - Oct. 11, 2013: Functions of ModelExchange and CoSimulation merged:
                      fmiInstantiateModel , fmiInstantiateSlave  -> fmiInstantiate
                      fmiFreeModelInstance, fmiFreeSlaveInstance -> fmiFreeInstance
                      fmiEnterModelInitializationMode, fmiEnterSlaveInitializationMode -> fmiEnterInitializationMode
                      fmiExitModelInitializationMode , fmiExitSlaveInitializationMode  -> fmiExitInitializationMode
                      fmiTerminateModel, fmiTerminateSlave  -> fmiTerminate
                      fmiResetSlave -> fmiReset (now also for ModelExchange and not only for CoSimulation)
                    Functions renamed:
                      fmiUpdateDiscreteStates -> fmiNewDiscreteStates
   - June 13, 2013: Functions removed:
                       fmiInitializeModel
                       fmiEventUpdate
                       fmiCompletedEventIteration
                       fmiInitializeSlave
                    Functions added:
                       fmiEnterModelInitializationMode
                       fmiExitModelInitializationMode
                       fmiEnterEventMode
                       fmiUpdateDiscreteStates
                       fmiEnterContinuousTimeMode
                       fmiEnterSlaveInitializationMode;
                       fmiExitSlaveInitializationMode;
   - Feb. 17, 2013: Portability improvements:
                       o DllExport changed to FMI_Export
                       o FUNCTION_PREFIX changed to FMI_FUNCTION_PREFIX
                       o Allow undefined FMI_FUNCTION_PREFIX (meaning no prefix is used)
                    Changed function name "fmiTerminate" to "fmiTerminateModel" (due to #113)
                    Changed function name "fmiGetNominalContinuousState" to
                                          "fmiGetNominalsOfContinuousStates"
                    Removed fmiGetStateValueReferences.
   - Nov. 14, 2011: Adapted to FMI 2.0:
                       o Split into two files (fmiFunctions.h, fmiTypes.h) in order
                         that code that dynamically loads an FMU can directly
                         utilize the header files).
                       o Added C++ encapsulation of C-part, in order that the header
                         file can be directly utilized in C++ code.
                       o fmiCallbackFunctions is passed as pointer to fmiInstantiateXXX
                       o stepFinished within fmiCallbackFunctions has as first
                         argument "fmiComponentEnvironment" and not "fmiComponent".
                       o New functions to get and set the complete FMU state
                         and to compute partial derivatives.
   - Nov.  4, 2010: Adapted to specification text:
                       o fmiGetModelTypesPlatform renamed to fmiGetTypesPlatform
                       o fmiInstantiateSlave: Argument GUID     replaced by fmuGUID
                                              Argument mimetype replaced by mimeType
                       o tabs replaced by spaces
   - Oct. 16, 2010: Functions for FMI for Co-simulation added
   - Jan. 20, 2010: stateValueReferencesChanged added to struct fmiEventInfo (ticket #27)
                    (by M. Otter, DLR)
                    Added WIN32 pragma to define the struct layout (ticket #34)
                    (by J. Mauss, QTronic)
   - Jan.  4, 2010: Removed argument intermediateResults from fmiInitialize
                    Renamed macro fmiGetModelFunctionsVersion to fmiGetVersion
                    Renamed macro fmiModelFunctionsVersion to fmiVersion
                    Replaced fmiModel by fmiComponent in decl of fmiInstantiateModel
                    (by J. Mauss, QTronic)
   - Dec. 17, 2009: Changed extension "me" to "fmi" (by Martin Otter, DLR).
   - Dez. 14, 2009: Added eventInfo to meInitialize and added
                    meGetNominalContinuousStates (by Martin Otter, DLR)
   - Sept. 9, 2009: Added DllExport (according to Peter Nilsson's suggestion)
                    (by A. Junghanns, QTronic)
   - Sept. 9, 2009: Changes according to FMI-meeting on July 21:
                    meInquireModelTypesVersion     -> meGetModelTypesPlatform
                    meInquireModelFunctionsVersion -> meGetModelFunctionsVersion
                    meSetStates                    -> meSetContinuousStates
                    meGetStates                    -> meGetContinuousStates
                    removal of meInitializeModelClass
                    removal of meGetTime
                    change of arguments of meInstantiateModel
                    change of arguments of meCompletedIntegratorStep
                    (by Martin Otter, DLR):
   - July 19, 2009: Added "me" as prefix to file names (by Martin Otter, DLR).
   - March 2, 2009: Changed function definitions according to the last design
                    meeting with additional improvements (by Martin Otter, DLR).
   - Dec. 3 , 2008: First version by Martin Otter (DLR) and Hans Olsson (Dynasim).


   Copyright (C) 2008-2011 MODELISAR consortium,
                 2012-2019 Modelica Association Project "FMI"
                 All rights reserved.

   This file is licensed by the copyright holders under the 2-Clause BSD License
   (https://opensource.org/licenses/BSD-2-Clause):

   ----------------------------------------------------------------------------
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

   - Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.

   - Redistributions in binary form must reproduce the above copyright notice,
     this list of conditions and the following disclaimer in the documentation
     and/or other materials provided with the distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
   OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
   WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
   OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
   ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
   ----------------------------------------------------------------------------
*/

#ifdef __cplusplus
extern "C" {
#endif

#ifndef fmi2TypesPlatform_h
#define fmi2TypesPlatform_h

/* Standard header file to define the argument types of the
   functions of the Functional Mock-up Interface 2.0.1
   This header file must be utilized both by the model and
   by the simulation engine.

   Revisions:
   - Sep. 29, 2019:  License changed to 2-clause BSD License (without extensions)
   - Apr.  9, 2014:  All prefixes "fmi" renamed to "fmi2" (decision from April 8)
   - Mar   31, 2014: New datatype fmiChar introduced.
   - Feb.  17, 2013: Changed fmiTypesPlatform from "standard32" to "default".
                     Removed fmiUndefinedValueReference since no longer needed
                     (because every state is defined in ScalarVariables).
   - March 20, 2012: Renamed from fmiPlatformTypes.h to fmiTypesPlatform.h
   - Nov.  14, 2011: Use the header file "fmiPlatformTypes.h" for FMI 2.0
                     both for "FMI for model exchange" and for "FMI for co-simulation"
                     New types "fmiComponentEnvironment", "fmiState", and "fmiByte".
                     The implementation of "fmiBoolean" is change from "char" to "int".
                     The #define "fmiPlatform" changed to "fmiTypesPlatform"
                     (in order that #define and function call are consistent)
   - Oct.   4, 2010: Renamed header file from "fmiModelTypes.h" to fmiPlatformTypes.h"
                     for the co-simulation interface
   - Jan.   4, 2010: Renamed meModelTypes_h to fmiModelTypes_h (by Mauss, QTronic)
   - Dec.  21, 2009: Changed "me" to "fmi" and "meModel" to "fmiComponent"
                     according to meeting on Dec. 18 (by Martin Otter, DLR)
   - Dec.   6, 2009: Added meUndefinedValueReference (by Martin Otter, DLR)
   - Sept.  9, 2009: Changes according to FMI-meeting on July 21:
                     Changed "version" to "platform", "standard" to "standard32",
                     Added a precise definition of "standard32" as comment
                     (by Martin Otter, DLR)
   - July  19, 2009: Added "me" as prefix to file names, added meTrue/meFalse,
                     and changed meValueReferenced from int to unsigned int
                     (by Martin Otter, DLR).
   - March  2, 2009: Moved enums and function pointer definitions to
                     ModelFunctions.h (by Martin Otter, DLR).
   - Dec.  3, 2008 : First version by Martin Otter (DLR) and
                     Hans Olsson (Dynasim).


   Copyright (C) 2008-2011 MODELISAR consortium,
                 2012-2019 Modelica Association Project "FMI"
                 All rights reserved.

   This file is licensed by the copyright holders under the 2-Clause BSD License
   (https://opensource.org/licenses/BSD-2-Clause):

   ----------------------------------------------------------------------------
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

   - Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.

   - Redistributions in binary form must reproduce the above copyright notice,
     this list of conditions and the following disclaimer in the documentation
     and/or other materials provided with the distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
   OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
   WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
   OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
   ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
   ----------------------------------------------------------------------------
*/

/* Platform (unique identification of this header file) */
#define fmi2TypesPlatform "default"

/* Type definitions of variables passed as arguments
   Version "default" means:

   fmi2Component           : an opaque object pointer
   fmi2ComponentEnvironment: an opaque object pointer
   fmi2FMUstate            : an opaque object pointer
   fmi2ValueReference      : handle to the value of a variable
   fmi2Real                : double precision floating-point data type
   fmi2Integer             : basic signed integer data type
   fmi2Boolean             : basic signed integer data type
   fmi2Char                : character data type
   fmi2String              : a pointer to a vector of fmi2Char characters
                             ('\0' terminated, UTF8 encoded)
   fmi2Byte                : smallest addressable unit of the machine, typically one byte.
*/
   typedef void*           fmi2Component;               /* Pointer to FMU instance       */
   typedef void*           fmi2ComponentEnvironment;    /* Pointer to FMU environment    */
   typedef void*           fmi2FMUstate;                /* Pointer to internal FMU state */
   typedef unsigned int    fmi2ValueReference;
   typedef double          fmi2Real   ;
   typedef int             fmi2Integer;
   typedef int             fmi2Boolean;
   typedef char            fmi2Char;
   typedef const fmi2Char* fmi2String;
   typedef char            fmi2Byte;

/* Values for fmi2Boolean  */
#define fmi2True  1
#define fmi2False 0


#endif /* fmi2TypesPlatform_h */
#ifndef fmi2FunctionTypes_h
#define fmi2FunctionTypes_h

#ifndef fmi2TypesPlatform_h
#define fmi2TypesPlatform_h

/* Standard header file to define the argument types of the
   functions of the Functional Mock-up Interface 2.0.1
   This header file must be utilized both by the model and
   by the simulation engine.

   Revisions:
   - Sep. 29, 2019:  License changed to 2-clause BSD License (without extensions)
   - Apr.  9, 2014:  All prefixes "fmi" renamed to "fmi2" (decision from April 8)
   - Mar   31, 2014: New datatype fmiChar introduced.
   - Feb.  17, 2013: Changed fmiTypesPlatform from "standard32" to "default".
                     Removed fmiUndefinedValueReference since no longer needed
                     (because every state is defined in ScalarVariables).
   - March 20, 2012: Renamed from fmiPlatformTypes.h to fmiTypesPlatform.h
   - Nov.  14, 2011: Use the header file "fmiPlatformTypes.h" for FMI 2.0
                     both for "FMI for model exchange" and for "FMI for co-simulation"
                     New types "fmiComponentEnvironment", "fmiState", and "fmiByte".
                     The implementation of "fmiBoolean" is change from "char" to "int".
                     The #define "fmiPlatform" changed to "fmiTypesPlatform"
                     (in order that #define and function call are consistent)
   - Oct.   4, 2010: Renamed header file from "fmiModelTypes.h" to fmiPlatformTypes.h"
                     for the co-simulation interface
   - Jan.   4, 2010: Renamed meModelTypes_h to fmiModelTypes_h (by Mauss, QTronic)
   - Dec.  21, 2009: Changed "me" to "fmi" and "meModel" to "fmiComponent"
                     according to meeting on Dec. 18 (by Martin Otter, DLR)
   - Dec.   6, 2009: Added meUndefinedValueReference (by Martin Otter, DLR)
   - Sept.  9, 2009: Changes according to FMI-meeting on July 21:
                     Changed "version" to "platform", "standard" to "standard32",
                     Added a precise definition of "standard32" as comment
                     (by Martin Otter, DLR)
   - July  19, 2009: Added "me" as prefix to file names, added meTrue/meFalse,
                     and changed meValueReferenced from int to unsigned int
                     (by Martin Otter, DLR).
   - March  2, 2009: Moved enums and function pointer definitions to
                     ModelFunctions.h (by Martin Otter, DLR).
   - Dec.  3, 2008 : First version by Martin Otter (DLR) and
                     Hans Olsson (Dynasim).


   Copyright (C) 2008-2011 MODELISAR consortium,
                 2012-2019 Modelica Association Project "FMI"
                 All rights reserved.

   This file is licensed by the copyright holders under the 2-Clause BSD License
   (https://opensource.org/licenses/BSD-2-Clause):

   ----------------------------------------------------------------------------
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

   - Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.

   - Redistributions in binary form must reproduce the above copyright notice,
     this list of conditions and the following disclaimer in the documentation
     and/or other materials provided with the distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
   OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
   WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
   OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
   ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
   ----------------------------------------------------------------------------
*/

/* Platform (unique identification of this header file) */
#define fmi2TypesPlatform "default"

/* Type definitions of variables passed as arguments
   Version "default" means:

   fmi2Component           : an opaque object pointer
   fmi2ComponentEnvironment: an opaque object pointer
   fmi2FMUstate            : an opaque object pointer
   fmi2ValueReference      : handle to the value of a variable
   fmi2Real                : double precision floating-point data type
   fmi2Integer             : basic signed integer data type
   fmi2Boolean             : basic signed integer data type
   fmi2Char                : character data type
   fmi2String              : a pointer to a vector of fmi2Char characters
                             ('\0' terminated, UTF8 encoded)
   fmi2Byte                : smallest addressable unit of the machine, typically one byte.
*/
   typedef void*           fmi2Component;               /* Pointer to FMU instance       */
   typedef void*           fmi2ComponentEnvironment;    /* Pointer to FMU environment    */
   typedef void*           fmi2FMUstate;                /* Pointer to internal FMU state */
   typedef unsigned int    fmi2ValueReference;
   typedef double          fmi2Real   ;
   typedef int             fmi2Integer;
   typedef int             fmi2Boolean;
   typedef char            fmi2Char;
   typedef const fmi2Char* fmi2String;
   typedef char            fmi2Byte;

/* Values for fmi2Boolean  */
#define fmi2True  1
#define fmi2False 0


#endif /* fmi2TypesPlatform_h */

/* This header file must be utilized when compiling an FMU or an FMI master.
   It declares data and function types for FMI 2.0.1

   Revisions:
   - Sep. 30, 2019: License changed to 2-clause BSD License (without extensions)
   - Jul.  5, 2019: Remove const modifier from fields of fmi2CallbackFunctions  (#216)
   - Sep.  6, 2018: Parameter names added to function prototypes
   - Apr.  9, 2014: all prefixes "fmi" renamed to "fmi2" (decision from April 8)
   - Apr.  3, 2014: Added #include <stddef.h> for size_t definition
   - Mar. 27, 2014: Added #include "fmiTypesPlatform.h" (#179)
   - Mar. 26, 2014: Introduced function argument "void" for the functions (#171)
                      fmiGetTypesPlatformTYPE and fmiGetVersionTYPE
   - Oct. 11, 2013: Functions of ModelExchange and CoSimulation merged:
                      fmiInstantiateModelTYPE , fmiInstantiateSlaveTYPE  -> fmiInstantiateTYPE
                      fmiFreeModelInstanceTYPE, fmiFreeSlaveInstanceTYPE -> fmiFreeInstanceTYPE
                      fmiEnterModelInitializationModeTYPE, fmiEnterSlaveInitializationModeTYPE -> fmiEnterInitializationModeTYPE
                      fmiExitModelInitializationModeTYPE , fmiExitSlaveInitializationModeTYPE  -> fmiExitInitializationModeTYPE
                      fmiTerminateModelTYPE , fmiTerminateSlaveTYPE  -> fmiTerminate
                      fmiResetSlave -> fmiReset (now also for ModelExchange and not only for CoSimulation)
                    Functions renamed
                      fmiUpdateDiscreteStatesTYPE -> fmiNewDiscreteStatesTYPE
                    Renamed elements of the enumeration fmiEventInfo
                      upcomingTimeEvent             -> nextEventTimeDefined // due to generic naming scheme: varDefined + var
                      newUpdateDiscreteStatesNeeded -> newDiscreteStatesNeeded;
   - June 13, 2013: Changed type fmiEventInfo
                    Functions removed:
                       fmiInitializeModelTYPE
                       fmiEventUpdateTYPE
                       fmiCompletedEventIterationTYPE
                       fmiInitializeSlaveTYPE
                    Functions added:
                       fmiEnterModelInitializationModeTYPE
                       fmiExitModelInitializationModeTYPE
                       fmiEnterEventModeTYPE
                       fmiUpdateDiscreteStatesTYPE
                       fmiEnterContinuousTimeModeTYPE
                       fmiEnterSlaveInitializationModeTYPE;
                       fmiExitSlaveInitializationModeTYPE;
   - Feb. 17, 2013: Added third argument to fmiCompletedIntegratorStepTYPE
                    Changed function name "fmiTerminateType" to "fmiTerminateModelType" (due to #113)
                    Changed function name "fmiGetNominalContinuousStateTYPE" to
                                          "fmiGetNominalsOfContinuousStatesTYPE"
                    Removed fmiGetStateValueReferencesTYPE.
   - Nov. 14, 2011: First public Version


   Copyright (C) 2008-2011 MODELISAR consortium,
                 2012-2019 Modelica Association Project "FMI"
                 All rights reserved.

   This file is licensed by the copyright holders under the 2-Clause BSD License
   (https://opensource.org/licenses/BSD-2-Clause):

   ----------------------------------------------------------------------------
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

   - Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.

   - Redistributions in binary form must reproduce the above copyright notice,
     this list of conditions and the following disclaimer in the documentation
     and/or other materials provided with the distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
   OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
   WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
   OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
   ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
   ----------------------------------------------------------------------------
*/

#ifdef __cplusplus
extern "C" {
#endif

/* make sure all compiler use the same alignment policies for structures */
#if defined _MSC_VER || defined __GNUC__
#pragma pack(push,8)
#endif

/* Include stddef.h, in order that size_t etc. is defined */
#include <stddef.h>


/* Type definitions */
typedef enum {
    fmi2OK,
    fmi2Warning,
    fmi2Discard,
    fmi2Error,
    fmi2Fatal,
    fmi2Pending
} fmi2Status;

typedef enum {
    fmi2ModelExchange,
    fmi2CoSimulation
} fmi2Type;

typedef enum {
    fmi2DoStepStatus,
    fmi2PendingStatus,
    fmi2LastSuccessfulTime,
    fmi2Terminated
} fmi2StatusKind;

typedef void      (*fmi2CallbackLogger)        (fmi2ComponentEnvironment componentEnvironment,
                                                fmi2String instanceName,
                                                fmi2Status status,
                                                fmi2String category,
                                                fmi2String message,
                                                ...);
typedef void*     (*fmi2CallbackAllocateMemory)(size_t nobj, size_t size);
typedef void      (*fmi2CallbackFreeMemory)    (void* obj);
typedef void      (*fmi2StepFinished)          (fmi2ComponentEnvironment componentEnvironment,
                                                fmi2Status status);

typedef struct {
   fmi2CallbackLogger         logger;
   fmi2CallbackAllocateMemory allocateMemory;
   fmi2CallbackFreeMemory     freeMemory;
   fmi2StepFinished           stepFinished;
   fmi2ComponentEnvironment   componentEnvironment;
} fmi2CallbackFunctions;

typedef struct {
   fmi2Boolean newDiscreteStatesNeeded;
   fmi2Boolean terminateSimulation;
   fmi2Boolean nominalsOfContinuousStatesChanged;
   fmi2Boolean valuesOfContinuousStatesChanged;
   fmi2Boolean nextEventTimeDefined;
   fmi2Real    nextEventTime;
} fmi2EventInfo;


/* reset alignment policy to the one set before reading this file */
#if defined _MSC_VER || defined __GNUC__
#pragma pack(pop)
#endif


/* Define fmi2 function pointer types to simplify dynamic loading */

/***************************************************
Types for Common Functions
****************************************************/

/* Inquire version numbers of header files and setting logging status */
    typedef const char* fmi2GetTypesPlatformTYPE(void);
    typedef const char* fmi2GetVersionTYPE(void);
    typedef fmi2Status  fmi2SetDebugLoggingTYPE(fmi2Component c,
                                                fmi2Boolean loggingOn,
                                                size_t nCategories,
                                                const fmi2String categories[]);

/* Creation and destruction of FMU instances and setting debug status */
    typedef fmi2Component fmi2InstantiateTYPE(fmi2String instanceName,
                                              fmi2Type fmuType,
                                              fmi2String fmuGUID,
                                              fmi2String fmuResourceLocation,
                                              const fmi2CallbackFunctions* functions,
                                              fmi2Boolean visible,
                                              fmi2Boolean loggingOn);
   typedef void          fmi2FreeInstanceTYPE(fmi2Component c);

/* Enter and exit initialization mode, terminate and reset */
    typedef fmi2Status fmi2SetupExperimentTYPE       (fmi2Component c,
                                                      fmi2Boolean toleranceDefined,
                                                      fmi2Real tolerance,
                                                      fmi2Real startTime,
                                                      fmi2Boolean stopTimeDefined,
                                                      fmi2Real stopTime);
   typedef fmi2Status fmi2EnterInitializationModeTYPE(fmi2Component c);
   typedef fmi2Status fmi2ExitInitializationModeTYPE (fmi2Component c);
   typedef fmi2Status fmi2TerminateTYPE              (fmi2Component c);
   typedef fmi2Status fmi2ResetTYPE                  (fmi2Component c);

/* Getting and setting variable values */
   typedef fmi2Status fmi2GetRealTYPE   (fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Real    value[]);
   typedef fmi2Status fmi2GetIntegerTYPE(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Integer value[]);
   typedef fmi2Status fmi2GetBooleanTYPE(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Boolean value[]);
   typedef fmi2Status fmi2GetStringTYPE (fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2String  value[]);

   typedef fmi2Status fmi2SetRealTYPE   (fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Real    value[]);
   typedef fmi2Status fmi2SetIntegerTYPE(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Integer value[]);
   typedef fmi2Status fmi2SetBooleanTYPE(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Boolean value[]);
   typedef fmi2Status fmi2SetStringTYPE (fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2String  value[]);

/* Getting and setting the internal FMU state */
   typedef fmi2Status fmi2GetFMUstateTYPE           (fmi2Component c, fmi2FMUstate* FMUstate);
   typedef fmi2Status fmi2SetFMUstateTYPE           (fmi2Component c, fmi2FMUstate  FMUstate);
   typedef fmi2Status fmi2FreeFMUstateTYPE          (fmi2Component c, fmi2FMUstate* FMUstate);
   typedef fmi2Status fmi2SerializedFMUstateSizeTYPE(fmi2Component c, fmi2FMUstate  FMUstate, size_t* size);
   typedef fmi2Status fmi2SerializeFMUstateTYPE     (fmi2Component c, fmi2FMUstate  FMUstate, fmi2Byte[], size_t size);
   typedef fmi2Status fmi2DeSerializeFMUstateTYPE   (fmi2Component c, const fmi2Byte serializedState[], size_t size, fmi2FMUstate* FMUstate);

/* Getting partial derivatives */
    typedef fmi2Status fmi2GetDirectionalDerivativeTYPE(fmi2Component c,
                                                        const fmi2ValueReference vUnknown_ref[], size_t nUnknown,
                                                        const fmi2ValueReference vKnown_ref[],   size_t nKnown,
                                                        const fmi2Real dvKnown[],
                                                        fmi2Real dvUnknown[]);

/***************************************************
Types for Functions for FMI2 for Model Exchange
****************************************************/

/* Enter and exit the different modes */
   typedef fmi2Status fmi2EnterEventModeTYPE         (fmi2Component c);
   typedef fmi2Status fmi2NewDiscreteStatesTYPE      (fmi2Component c, fmi2EventInfo* fmi2eventInfo);
   typedef fmi2Status fmi2EnterContinuousTimeModeTYPE(fmi2Component c);
   typedef fmi2Status fmi2CompletedIntegratorStepTYPE(fmi2Component c,
                                                      fmi2Boolean   noSetFMUStatePriorToCurrentPoint,
                                                      fmi2Boolean*  enterEventMode,
                                                      fmi2Boolean*  terminateSimulation);

/* Providing independent variables and re-initialization of caching */
   typedef fmi2Status fmi2SetTimeTYPE            (fmi2Component c, fmi2Real time);
   typedef fmi2Status fmi2SetContinuousStatesTYPE(fmi2Component c, const fmi2Real x[], size_t nx);

/* Evaluation of the model equations */
   typedef fmi2Status fmi2GetDerivativesTYPE               (fmi2Component c, fmi2Real derivatives[],     size_t nx);
   typedef fmi2Status fmi2GetEventIndicatorsTYPE           (fmi2Component c, fmi2Real eventIndicators[], size_t ni);
   typedef fmi2Status fmi2GetContinuousStatesTYPE          (fmi2Component c, fmi2Real x[],               size_t nx);
   typedef fmi2Status fmi2GetNominalsOfContinuousStatesTYPE(fmi2Component c, fmi2Real x_nominal[],       size_t nx);


/***************************************************
Types for Functions for FMI2 for Co-Simulation
****************************************************/

/* Simulating the slave */
    typedef fmi2Status fmi2SetRealInputDerivativesTYPE (fmi2Component c,
                                                        const fmi2ValueReference vr[], size_t nvr,
                                                        const fmi2Integer order[],
                                                        const fmi2Real value[]);
    typedef fmi2Status fmi2GetRealOutputDerivativesTYPE(fmi2Component c,
                                                        const fmi2ValueReference vr[], size_t nvr,
                                                        const fmi2Integer order[],
                                                        fmi2Real value[]);
    typedef fmi2Status fmi2DoStepTYPE   (fmi2Component c,
                                         fmi2Real      currentCommunicationPoint,
                                         fmi2Real      communicationStepSize,
                                         fmi2Boolean   noSetFMUStatePriorToCurrentPoint);
   typedef fmi2Status fmi2CancelStepTYPE(fmi2Component c);

/* Inquire slave status */
   typedef fmi2Status fmi2GetStatusTYPE       (fmi2Component c, const fmi2StatusKind s, fmi2Status*  value);
   typedef fmi2Status fmi2GetRealStatusTYPE   (fmi2Component c, const fmi2StatusKind s, fmi2Real*    value);
   typedef fmi2Status fmi2GetIntegerStatusTYPE(fmi2Component c, const fmi2StatusKind s, fmi2Integer* value);
   typedef fmi2Status fmi2GetBooleanStatusTYPE(fmi2Component c, const fmi2StatusKind s, fmi2Boolean* value);
   typedef fmi2Status fmi2GetStringStatusTYPE (fmi2Component c, const fmi2StatusKind s, fmi2String*  value);


#ifdef __cplusplus
}  /* end of extern "C" { */
#endif

#endif /* fmi2FunctionTypes_h */
#include <stdlib.h>


/*
  Export FMI2 API functions on Windows and under GCC.
  If custom linking is desired then the FMI2_Export must be
  defined before including this file. For instance,
  it may be set to __declspec(dllimport).
*/
#if !defined(FMI2_Export)
  #if !defined(FMI2_FUNCTION_PREFIX)
    #if defined _WIN32 || defined __CYGWIN__
     /* Note: both gcc & MSVC on Windows support this syntax. */
        #define FMI2_Export __declspec(dllexport)
    #else
      #if __GNUC__ >= 4
        #define FMI2_Export __attribute__ ((visibility ("default")))
      #else
        #define FMI2_Export
      #endif
    #endif
  #else
    #define FMI2_Export
  #endif
#endif

/* Macros to construct the real function name
   (prepend function name by FMI2_FUNCTION_PREFIX) */
#if defined(FMI2_FUNCTION_PREFIX)
  #define fmi2Paste(a,b)     a ## b
  #define fmi2PasteB(a,b)    fmi2Paste(a,b)
  #define fmi2FullName(name) fmi2PasteB(FMI2_FUNCTION_PREFIX, name)
#else
  #define fmi2FullName(name) name
#endif

/***************************************************
Common Functions
****************************************************/
#define fmi2GetTypesPlatform         fmi2FullName(fmi2GetTypesPlatform)
#define fmi2GetVersion               fmi2FullName(fmi2GetVersion)
#define fmi2SetDebugLogging          fmi2FullName(fmi2SetDebugLogging)
#define fmi2Instantiate              fmi2FullName(fmi2Instantiate)
#define fmi2FreeInstance             fmi2FullName(fmi2FreeInstance)
#define fmi2SetupExperiment          fmi2FullName(fmi2SetupExperiment)
#define fmi2EnterInitializationMode  fmi2FullName(fmi2EnterInitializationMode)
#define fmi2ExitInitializationMode   fmi2FullName(fmi2ExitInitializationMode)
#define fmi2Terminate                fmi2FullName(fmi2Terminate)
#define fmi2Reset                    fmi2FullName(fmi2Reset)
#define fmi2GetReal                  fmi2FullName(fmi2GetReal)
#define fmi2GetInteger               fmi2FullName(fmi2GetInteger)
#define fmi2GetBoolean               fmi2FullName(fmi2GetBoolean)
#define fmi2GetString                fmi2FullName(fmi2GetString)
#define fmi2SetReal                  fmi2FullName(fmi2SetReal)
#define fmi2SetInteger               fmi2FullName(fmi2SetInteger)
#define fmi2SetBoolean               fmi2FullName(fmi2SetBoolean)
#define fmi2SetString                fmi2FullName(fmi2SetString)
#define fmi2GetFMUstate              fmi2FullName(fmi2GetFMUstate)
#define fmi2SetFMUstate              fmi2FullName(fmi2SetFMUstate)
#define fmi2FreeFMUstate             fmi2FullName(fmi2FreeFMUstate)
#define fmi2SerializedFMUstateSize   fmi2FullName(fmi2SerializedFMUstateSize)
#define fmi2SerializeFMUstate        fmi2FullName(fmi2SerializeFMUstate)
#define fmi2DeSerializeFMUstate      fmi2FullName(fmi2DeSerializeFMUstate)
#define fmi2GetDirectionalDerivative fmi2FullName(fmi2GetDirectionalDerivative)


/***************************************************
Functions for FMI2 for Model Exchange
****************************************************/
#define fmi2EnterEventMode                fmi2FullName(fmi2EnterEventMode)
#define fmi2NewDiscreteStates             fmi2FullName(fmi2NewDiscreteStates)
#define fmi2EnterContinuousTimeMode       fmi2FullName(fmi2EnterContinuousTimeMode)
#define fmi2CompletedIntegratorStep       fmi2FullName(fmi2CompletedIntegratorStep)
#define fmi2SetTime                       fmi2FullName(fmi2SetTime)
#define fmi2SetContinuousStates           fmi2FullName(fmi2SetContinuousStates)
#define fmi2GetDerivatives                fmi2FullName(fmi2GetDerivatives)
#define fmi2GetEventIndicators            fmi2FullName(fmi2GetEventIndicators)
#define fmi2GetContinuousStates           fmi2FullName(fmi2GetContinuousStates)
#define fmi2GetNominalsOfContinuousStates fmi2FullName(fmi2GetNominalsOfContinuousStates)


/***************************************************
Functions for FMI2 for Co-Simulation
****************************************************/
#define fmi2SetRealInputDerivatives      fmi2FullName(fmi2SetRealInputDerivatives)
#define fmi2GetRealOutputDerivatives     fmi2FullName(fmi2GetRealOutputDerivatives)
#define fmi2DoStep                       fmi2FullName(fmi2DoStep)
#define fmi2CancelStep                   fmi2FullName(fmi2CancelStep)
#define fmi2GetStatus                    fmi2FullName(fmi2GetStatus)
#define fmi2GetRealStatus                fmi2FullName(fmi2GetRealStatus)
#define fmi2GetIntegerStatus             fmi2FullName(fmi2GetIntegerStatus)
#define fmi2GetBooleanStatus             fmi2FullName(fmi2GetBooleanStatus)
#define fmi2GetStringStatus              fmi2FullName(fmi2GetStringStatus)

/* Version number */
#define fmi2Version "2.0"


/***************************************************
Common Functions
****************************************************/

/* Inquire version numbers of header files */
   FMI2_Export fmi2GetTypesPlatformTYPE fmi2GetTypesPlatform;
   FMI2_Export fmi2GetVersionTYPE       fmi2GetVersion;
   FMI2_Export fmi2SetDebugLoggingTYPE  fmi2SetDebugLogging;

/* Creation and destruction of FMU instances */
   FMI2_Export fmi2InstantiateTYPE  fmi2Instantiate;
   FMI2_Export fmi2FreeInstanceTYPE fmi2FreeInstance;

/* Enter and exit initialization mode, terminate and reset */
   FMI2_Export fmi2SetupExperimentTYPE         fmi2SetupExperiment;
   FMI2_Export fmi2EnterInitializationModeTYPE fmi2EnterInitializationMode;
   FMI2_Export fmi2ExitInitializationModeTYPE  fmi2ExitInitializationMode;
   FMI2_Export fmi2TerminateTYPE               fmi2Terminate;
   FMI2_Export fmi2ResetTYPE                   fmi2Reset;

/* Getting and setting variables values */
   FMI2_Export fmi2GetRealTYPE    fmi2GetReal;
   FMI2_Export fmi2GetIntegerTYPE fmi2GetInteger;
   FMI2_Export fmi2GetBooleanTYPE fmi2GetBoolean;
   FMI2_Export fmi2GetStringTYPE  fmi2GetString;

   FMI2_Export fmi2SetRealTYPE    fmi2SetReal;
   FMI2_Export fmi2SetIntegerTYPE fmi2SetInteger;
   FMI2_Export fmi2SetBooleanTYPE fmi2SetBoolean;
   FMI2_Export fmi2SetStringTYPE  fmi2SetString;

/* Getting and setting the internal FMU state */
   FMI2_Export fmi2GetFMUstateTYPE            fmi2GetFMUstate;
   FMI2_Export fmi2SetFMUstateTYPE            fmi2SetFMUstate;
   FMI2_Export fmi2FreeFMUstateTYPE           fmi2FreeFMUstate;
   FMI2_Export fmi2SerializedFMUstateSizeTYPE fmi2SerializedFMUstateSize;
   FMI2_Export fmi2SerializeFMUstateTYPE      fmi2SerializeFMUstate;
   FMI2_Export fmi2DeSerializeFMUstateTYPE    fmi2DeSerializeFMUstate;

/* Getting partial derivatives */
   FMI2_Export fmi2GetDirectionalDerivativeTYPE fmi2GetDirectionalDerivative;


/***************************************************
Functions for FMI2 for Model Exchange
****************************************************/

/* Enter and exit the different modes */
   FMI2_Export fmi2EnterEventModeTYPE               fmi2EnterEventMode;
   FMI2_Export fmi2NewDiscreteStatesTYPE            fmi2NewDiscreteStates;
   FMI2_Export fmi2EnterContinuousTimeModeTYPE      fmi2EnterContinuousTimeMode;
   FMI2_Export fmi2CompletedIntegratorStepTYPE      fmi2CompletedIntegratorStep;

/* Providing independent variables and re-initialization of caching */
   FMI2_Export fmi2SetTimeTYPE             fmi2SetTime;
   FMI2_Export fmi2SetContinuousStatesTYPE fmi2SetContinuousStates;

/* Evaluation of the model equations */
   FMI2_Export fmi2GetDerivativesTYPE                fmi2GetDerivatives;
   FMI2_Export fmi2GetEventIndicatorsTYPE            fmi2GetEventIndicators;
   FMI2_Export fmi2GetContinuousStatesTYPE           fmi2GetContinuousStates;
   FMI2_Export fmi2GetNominalsOfContinuousStatesTYPE fmi2GetNominalsOfContinuousStates;


/***************************************************
Functions for FMI2 for Co-Simulation
****************************************************/

/* Simulating the slave */
   FMI2_Export fmi2SetRealInputDerivativesTYPE  fmi2SetRealInputDerivatives;
   FMI2_Export fmi2GetRealOutputDerivativesTYPE fmi2GetRealOutputDerivatives;

   FMI2_Export fmi2DoStepTYPE     fmi2DoStep;
   FMI2_Export fmi2CancelStepTYPE fmi2CancelStep;

/* Inquire slave status */
   FMI2_Export fmi2GetStatusTYPE        fmi2GetStatus;
   FMI2_Export fmi2GetRealStatusTYPE    fmi2GetRealStatus;
   FMI2_Export fmi2GetIntegerStatusTYPE fmi2GetIntegerStatus;
   FMI2_Export fmi2GetBooleanStatusTYPE fmi2GetBooleanStatus;
   FMI2_Export fmi2GetStringStatusTYPE  fmi2GetStringStatus;

#ifdef __cplusplus
}  /* end of extern "C" { */
#endif

#endif /* fmi2Functions_h */
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <Windows.h>
#endif

#include <stdbool.h>
#include <stddef.h>

#ifndef FMI_MAX_MESSAGE_LENGTH
#define FMI_MAX_MESSAGE_LENGTH 4096
#endif

#ifndef FMI_STATIC
#define FMI_STATIC
#endif

typedef enum {
    FMIOK,
    FMIWarning,
    FMIDiscard,
    FMIError,
    FMIFatal,
    FMIPending
} FMIStatus;

typedef enum {

    // FMI 3.0 variable types
    FMIFloat32Type,
    FMIDiscreteFloat32Type,
    FMIFloat64Type,
    FMIDiscreteFloat64Type,
    FMIInt8Type,
    FMIUInt8Type,
    FMIInt16Type,
    FMIUInt16Type,
    FMIInt32Type,
    FMIUInt32Type,
    FMIInt64Type,
    FMIUInt64Type,
    FMIBooleanType,
    FMIStringType,
    FMIBinaryType,
    FMIClockType,

    // non-variable types used as parameters
    FMIValueReferenceType,
    FMISizeTType,

    // aliases for FMI 1.0 and 2.0 variable types
    FMIRealType         = FMIFloat64Type,
    FMIDiscreteRealType = FMIDiscreteFloat64Type,
    FMIIntegerType      = FMIInt32Type,

} FMIVariableType;

typedef enum {
    FMIVersion1 = 1,
    FMIVersion2 = 2,
    FMIVersion3 = 3
} FMIVersion;

typedef enum {
    FMIModelExchange,
    FMICoSimulation,
    FMIScheduledExecution
} FMIInterfaceType;

typedef enum {
    FMI2StartAndEndState        = 1 << 0,
    FMI2InstantiatedState       = 1 << 1,
    FMI2InitializationModeState = 1 << 2,

    // model exchange states
    FMI2EventModeState          = 1 << 3,
    FMI2ContinuousTimeModeState = 1 << 4,

    // co-simulation states
    FMI2StepCompleteState       = 1 << 5,
    FMI2StepInProgressState     = 1 << 6,
    FMI2StepFailedState         = 1 << 7,
    FMI2StepCanceledState       = 1 << 8,

    FMI2TerminatedState         = 1 << 9,
    FMI2ErrorState              = 1 << 10,
    FMI2FatalState              = 1 << 11,
} FMI2State;

typedef unsigned int FMIValueReference;

typedef struct FMIInstance_ FMIInstance;

typedef struct FMI1Functions_ FMI1Functions;

typedef struct FMI2Functions_ FMI2Functions;

typedef struct FMI3Functions_ FMI3Functions;

typedef void FMILogFunctionCall(FMIInstance *instance, FMIStatus status, const char *message);

typedef void FMILogMessage(FMIInstance *instance, FMIStatus status, const char *category, const char *message);

struct FMIInstance_ {

    FMI1Functions *fmi1Functions;
    FMI2Functions *fmi2Functions;
    FMI3Functions *fmi3Functions;

#ifdef _WIN32
    HMODULE libraryHandle;
#else
    void *libraryHandle;
#endif

    void *userData;

    FMILogMessage      *logMessage;
    FMILogFunctionCall *logFunctionCall;

    double time;

    char* logMessageBuffer;
    size_t logMessageBufferSize;
    size_t logMessageBufferPosition;

    void *component;

    const char *name;

    bool logFMICalls;

    FMI2State state;

    FMIStatus status;

    FMIVersion fmiVersion;

    FMIInterfaceType interfaceType;

};

FMI_STATIC FMIInstance* FMICreateInstance(const char* instanceName, const char* libraryPath, FMILogMessage* logMessage, FMILogFunctionCall* logFunctionCall);

FMI_STATIC void FMIFreeInstance(FMIInstance *instance);

FMI_STATIC void FMIClearLogMessageBuffer(FMIInstance* instance);

FMI_STATIC void FMIAppendToLogMessageBuffer(FMIInstance* instance, const char* format, ...);

FMI_STATIC void FMIAppendArrayToLogMessageBuffer(FMIInstance* instance, const void* values, size_t nValues, const size_t sizes[], FMIVariableType variableType);

FMI_STATIC FMIStatus FMIURIToPath(const char *uri, char *path, const size_t pathLength);

FMI_STATIC FMIStatus FMIPathToURI(const char *path, char *uri, const size_t uriLength);

FMI_STATIC FMIStatus FMIPlatformBinaryPath(const char *unzipdir, const char *modelIdentifier, FMIVersion fmiVersion, char *platformBinaryPath, size_t size);

#ifdef __cplusplus
}  /* end of extern "C" { */
#endif

struct FMI2Functions_ {

    fmi2CallbackFunctions            callbacks;
    fmi2EventInfo                    eventInfo;

    /***************************************************
    Common Functions for FMI 2.0
    ****************************************************/

    /* required functions */
    fmi2GetTypesPlatformTYPE         *fmi2GetTypesPlatform;
    fmi2GetVersionTYPE               *fmi2GetVersion;
    fmi2SetDebugLoggingTYPE          *fmi2SetDebugLogging;
    fmi2InstantiateTYPE              *fmi2Instantiate;
    fmi2FreeInstanceTYPE             *fmi2FreeInstance;
    fmi2SetupExperimentTYPE          *fmi2SetupExperiment;
    fmi2EnterInitializationModeTYPE  *fmi2EnterInitializationMode;
    fmi2ExitInitializationModeTYPE   *fmi2ExitInitializationMode;
    fmi2TerminateTYPE                *fmi2Terminate;
    fmi2ResetTYPE                    *fmi2Reset;
    fmi2GetRealTYPE                  *fmi2GetReal;
    fmi2GetIntegerTYPE               *fmi2GetInteger;
    fmi2GetBooleanTYPE               *fmi2GetBoolean;
    fmi2GetStringTYPE                *fmi2GetString;
    fmi2SetRealTYPE                  *fmi2SetReal;
    fmi2SetIntegerTYPE               *fmi2SetInteger;
    fmi2SetBooleanTYPE               *fmi2SetBoolean;
    fmi2SetStringTYPE                *fmi2SetString;

    /* optional functions */
    fmi2GetFMUstateTYPE              *fmi2GetFMUstate;
    fmi2SetFMUstateTYPE              *fmi2SetFMUstate;
    fmi2FreeFMUstateTYPE             *fmi2FreeFMUstate;
    fmi2SerializedFMUstateSizeTYPE   *fmi2SerializedFMUstateSize;
    fmi2SerializeFMUstateTYPE        *fmi2SerializeFMUstate;
    fmi2DeSerializeFMUstateTYPE      *fmi2DeSerializeFMUstate;
    fmi2GetDirectionalDerivativeTYPE *fmi2GetDirectionalDerivative;

    /***************************************************
    Functions for FMI 2.0 for Model Exchange
    ****************************************************/

    fmi2EnterEventModeTYPE                *fmi2EnterEventMode;
    fmi2NewDiscreteStatesTYPE             *fmi2NewDiscreteStates;
    fmi2EnterContinuousTimeModeTYPE       *fmi2EnterContinuousTimeMode;
    fmi2CompletedIntegratorStepTYPE       *fmi2CompletedIntegratorStep;
    fmi2SetTimeTYPE                       *fmi2SetTime;
    fmi2SetContinuousStatesTYPE           *fmi2SetContinuousStates;
    fmi2GetDerivativesTYPE                *fmi2GetDerivatives;
    fmi2GetEventIndicatorsTYPE            *fmi2GetEventIndicators;
    fmi2GetContinuousStatesTYPE           *fmi2GetContinuousStates;
    fmi2GetNominalsOfContinuousStatesTYPE *fmi2GetNominalsOfContinuousStates;

    /***************************************************
    Functions for FMI 2.0 for Co-Simulation
    ****************************************************/

    fmi2SetRealInputDerivativesTYPE  *fmi2SetRealInputDerivatives;
    fmi2GetRealOutputDerivativesTYPE *fmi2GetRealOutputDerivatives;
    fmi2DoStepTYPE                   *fmi2DoStep;
    fmi2CancelStepTYPE               *fmi2CancelStep;
    fmi2GetStatusTYPE                *fmi2GetStatus;
    fmi2GetRealStatusTYPE            *fmi2GetRealStatus;
    fmi2GetIntegerStatusTYPE         *fmi2GetIntegerStatus;
    fmi2GetBooleanStatusTYPE         *fmi2GetBooleanStatus;
    fmi2GetStringStatusTYPE          *fmi2GetStringStatus;

};


/***************************************************
Common Functions
****************************************************/

/* Inquire version numbers of header files and setting logging status */
FMI_STATIC const char* FMI2GetTypesPlatform(FMIInstance *instance);

FMI_STATIC const char* FMI2GetVersion(FMIInstance *instance);

FMI_STATIC FMIStatus FMI2SetDebugLogging(FMIInstance *instance, fmi2Boolean loggingOn, size_t nCategories, const fmi2String categories[]);

FMI_STATIC FMIStatus FMI2Instantiate(FMIInstance *instance, const char *fmuResourceLocation, fmi2Type fmuType, fmi2String fmuGUID,
    fmi2Boolean visible, fmi2Boolean loggingOn);

FMI_STATIC void FMI2FreeInstance(FMIInstance *instance);

/* Enter and exit initialization mode, terminate and reset */
FMI_STATIC FMIStatus FMI2SetupExperiment(FMIInstance *instance,
    fmi2Boolean toleranceDefined,
    fmi2Real tolerance,
    fmi2Real startTime,
    fmi2Boolean stopTimeDefined,
    fmi2Real stopTime);

FMI_STATIC FMIStatus FMI2EnterInitializationMode(FMIInstance *instance);

FMI_STATIC FMIStatus FMI2ExitInitializationMode(FMIInstance *instance);

FMI_STATIC FMIStatus FMI2Terminate(FMIInstance *instance);

FMI_STATIC FMIStatus FMI2Reset(FMIInstance *instance);

FMI_STATIC FMIStatus FMI2SetupExperiment(FMIInstance *instance,
    fmi2Boolean toleranceDefined,
    fmi2Real tolerance,
    fmi2Real startTime,
    fmi2Boolean stopTimeDefined,
    fmi2Real stopTime);

/* Getting and setting variable values */
FMI_STATIC FMIStatus FMI2GetReal(FMIInstance *instance, const fmi2ValueReference vr[], size_t nvr, fmi2Real value[]);

FMI_STATIC FMIStatus FMI2GetInteger(FMIInstance *instance, const fmi2ValueReference vr[], size_t nvr, fmi2Integer value[]);

FMI_STATIC FMIStatus FMI2GetBoolean(FMIInstance *instance, const fmi2ValueReference vr[], size_t nvr, fmi2Boolean value[]);

FMI_STATIC FMIStatus FMI2GetString(FMIInstance *instance, const fmi2ValueReference vr[], size_t nvr, fmi2String  value[]);

FMI_STATIC FMIStatus FMI2SetReal(FMIInstance *instance, const fmi2ValueReference vr[], size_t nvr, const fmi2Real    value[]);

FMI_STATIC FMIStatus FMI2SetInteger(FMIInstance *instance, const fmi2ValueReference vr[], size_t nvr, const fmi2Integer value[]);

FMI_STATIC FMIStatus FMI2SetBoolean(FMIInstance *instance, const fmi2ValueReference vr[], size_t nvr, const fmi2Boolean value[]);

FMI_STATIC FMIStatus FMI2SetString(FMIInstance *instance, const fmi2ValueReference vr[], size_t nvr, const fmi2String  value[]);

/* Getting and setting the internal FMU state */
FMI_STATIC FMIStatus FMI2GetFMUstate(FMIInstance *instance, fmi2FMUstate* FMUstate);

FMI_STATIC FMIStatus FMI2SetFMUstate(FMIInstance *instance, fmi2FMUstate  FMUstate);

FMI_STATIC FMIStatus FMI2FreeFMUstate(FMIInstance *instance, fmi2FMUstate* FMUstate);

FMI_STATIC FMIStatus FMI2SerializedFMUstateSize(FMIInstance *instance, fmi2FMUstate  FMUstate, size_t* size);

FMI_STATIC FMIStatus FMI2SerializeFMUstate(FMIInstance *instance, fmi2FMUstate  FMUstate, fmi2Byte serializedState[], size_t size);

FMI_STATIC FMIStatus FMI2DeSerializeFMUstate(FMIInstance *instance, const fmi2Byte serializedState[], size_t size, fmi2FMUstate* FMUstate);

/* Getting partial derivatives */
FMI_STATIC FMIStatus FMI2GetDirectionalDerivative(FMIInstance *instance,
    const fmi2ValueReference vUnknown_ref[], size_t nUnknown,
    const fmi2ValueReference vKnown_ref[], size_t nKnown,
    const fmi2Real dvKnown[],
    fmi2Real dvUnknown[]);

/***************************************************
Types for Functions for FMI2 for Model Exchange
****************************************************/

/* Enter and exit the different modes */
FMI_STATIC FMIStatus FMI2EnterEventMode(FMIInstance *instance);

FMI_STATIC FMIStatus FMI2NewDiscreteStates(FMIInstance *instance, fmi2EventInfo *eventInfo);

FMI_STATIC FMIStatus FMI2EnterContinuousTimeMode(FMIInstance *instance);

FMI_STATIC FMIStatus FMI2CompletedIntegratorStep(FMIInstance *instance,
    fmi2Boolean   noSetFMUStatePriorToCurrentPoint,
    fmi2Boolean*  enterEventMode,
    fmi2Boolean*  terminateSimulation);

/* Providing independent variables and re-initialization of caching */
FMI_STATIC FMIStatus FMI2SetTime(FMIInstance *instance, fmi2Real time);

FMI_STATIC FMIStatus FMI2SetContinuousStates(FMIInstance *instance, const fmi2Real x[], size_t nx);

/* Evaluation of the model equations */
FMI_STATIC FMIStatus FMI2GetDerivatives(FMIInstance *instance, fmi2Real derivatives[], size_t nx);

FMI_STATIC FMIStatus FMI2GetEventIndicators(FMIInstance *instance, fmi2Real eventIndicators[], size_t ni);

FMI_STATIC FMIStatus FMI2GetContinuousStates(FMIInstance *instance, fmi2Real x[], size_t nx);

FMI_STATIC FMIStatus FMI2GetNominalsOfContinuousStates(FMIInstance *instance, fmi2Real x_nominal[], size_t nx);


/***************************************************
Types for Functions for FMI2 for Co-Simulation
****************************************************/

/* Simulating the slave */
FMI_STATIC FMIStatus FMI2SetRealInputDerivatives(FMIInstance *instance,
    const fmi2ValueReference vr[], size_t nvr,
    const fmi2Integer order[],
    const fmi2Real value[]);

FMI_STATIC FMIStatus FMI2GetRealOutputDerivatives(FMIInstance *instance,
    const fmi2ValueReference vr[], size_t nvr,
    const fmi2Integer order[],
    fmi2Real value[]);

FMI_STATIC FMIStatus FMI2DoStep(FMIInstance *instance,
    fmi2Real      currentCommunicationPoint,
    fmi2Real      communicationStepSize,
    fmi2Boolean   noSetFMUStatePriorToCurrentPoint);

FMI_STATIC FMIStatus FMI2CancelStep(FMIInstance *instance);

/* Inquire slave status */
FMI_STATIC FMIStatus FMI2GetStatus(FMIInstance *instance, const fmi2StatusKind s, fmi2Status* value);

FMI_STATIC FMIStatus FMI2GetRealStatus(FMIInstance *instance, const fmi2StatusKind s, fmi2Real* value);

FMI_STATIC FMIStatus FMI2GetIntegerStatus(FMIInstance *instance, const fmi2StatusKind s, fmi2Integer* value);

FMI_STATIC FMIStatus FMI2GetBooleanStatus(FMIInstance *instance, const fmi2StatusKind s, fmi2Boolean* value);

FMI_STATIC FMIStatus FMI2GetStringStatus(FMIInstance *instance, const fmi2StatusKind s, fmi2String*  value);

#ifdef __cplusplus
}  /* end of extern "C" { */
#endif


static void cb_logMessage2(fmi2ComponentEnvironment componentEnvironment, fmi2String instanceName, fmi2Status status, fmi2String category, fmi2String message, ...) {

    if (!componentEnvironment) return;

    FMIInstance *instance = componentEnvironment;

    char buf[FMI_MAX_MESSAGE_LENGTH];

    va_list args;

    va_start(args, message);
    vsnprintf(buf, FMI_MAX_MESSAGE_LENGTH, message, args);
    va_end(args);

    if (!instance->logMessage) return;

    instance->logMessage(instance, (FMIStatus)status, category, buf);
}

#if defined(FMI2_FUNCTION_PREFIX)
#define LOAD_SYMBOL(f) \
do { \
    instance->fmi2Functions->fmi2 ## f = fmi2 ## f; \
} while (0)
#elif defined(_WIN32)
#define LOAD_SYMBOL(f) \
do { \
    instance->fmi2Functions->fmi2 ## f = (fmi2 ## f ## TYPE*)GetProcAddress(instance->libraryHandle, "fmi2" #f); \
    if (!instance->fmi2Functions->fmi2 ## f) { \
        instance->logMessage(instance, FMIFatal, "fatal", "Symbol fmi2" #f " is missing in shared library."); \
        return fmi2Fatal; \
    }\
} while (0)
#else
#define LOAD_SYMBOL(f) \
do { \
    instance->fmi2Functions->fmi2 ## f = (fmi2 ## f ## TYPE*)dlsym(instance->libraryHandle, "fmi2" #f); \
    if (!instance->fmi2Functions->fmi2 ## f) { \
        instance->logMessage(instance, FMIFatal, "fatal", "Symbol fmi2" #f " is missing in shared library."); \
        return FMIFatal; \
    } \
} while (0)
#endif

#define CALL(f) \
do { \
    const FMIStatus status = (FMIStatus)instance->fmi2Functions->fmi2 ## f (instance->component); \
    if (instance->logFunctionCall) { \
        instance->logFunctionCall(instance, status, "fmi2" #f "()"); \
    } \
    return status; \
} while (0)

#define CALL_ARGS(f, m, ...) \
do { \
    const FMIStatus status = (FMIStatus)instance->fmi2Functions-> fmi2 ## f (instance->component, __VA_ARGS__); \
    if (instance->logFunctionCall) { \
        FMIClearLogMessageBuffer(instance); \
        FMIAppendToLogMessageBuffer(instance, "fmi2" #f "(" m ")", __VA_ARGS__); \
        instance->logFunctionCall(instance, status, instance->logMessageBuffer); \
    } \
    return status; \
} while (0)

#define CALL_ARRAY(s, t) \
do { \
    const FMIStatus status = (FMIStatus)instance->fmi2Functions->fmi2 ## s ## t(instance->component, vr, nvr, value); \
    if (instance->logFunctionCall) { \
        FMIClearLogMessageBuffer(instance); \
        FMIAppendToLogMessageBuffer(instance, "fmi2" #s #t "(vr={"); \
        FMIAppendArrayToLogMessageBuffer(instance, vr, nvr, NULL, FMIValueReferenceType); \
        FMIAppendToLogMessageBuffer(instance, "}, nvr=%zu, value={", nvr); \
        FMIAppendArrayToLogMessageBuffer(instance, value, nvr, NULL, FMI ## t ## Type); \
        FMIAppendToLogMessageBuffer(instance, "})"); \
        instance->logFunctionCall(instance, status, instance->logMessageBuffer); \
    } \
    return status; \
} while (0)

/***************************************************
Common Functions
****************************************************/

/* Inquire version numbers of header files and setting logging status */
const char* FMI2GetTypesPlatform(FMIInstance *instance) {

    if (instance->logFunctionCall) {
        instance->logFunctionCall(instance, FMIOK, "fmi2GetTypesPlatform()");
    }

    return instance->fmi2Functions->fmi2GetTypesPlatform();
}

const char* FMI2GetVersion(FMIInstance *instance) {

    if (instance->logFunctionCall) {
        instance->logFunctionCall(instance, FMIOK, "fmi2GetVersion()");
    }

    return instance->fmi2Functions->fmi2GetVersion();
}

FMIStatus FMI2SetDebugLogging(FMIInstance *instance, fmi2Boolean loggingOn, size_t nCategories, const fmi2String categories[]) {

    const FMIStatus status = (FMIStatus)instance->fmi2Functions->fmi2SetDebugLogging(instance->component, loggingOn, nCategories, categories);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi2SetDebugLogging(loggingOn=%d, nCategories=%zu, categories={");
        FMIAppendArrayToLogMessageBuffer(instance, categories, nCategories, NULL, FMIStringType);
        FMIAppendToLogMessageBuffer(instance, "})");
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

/* Creation and destruction of FMU instances and setting debug status */
FMIStatus FMI2Instantiate(FMIInstance *instance, const char *fmuResourceLocation, fmi2Type fmuType, fmi2String fmuGUID,
    fmi2Boolean visible, fmi2Boolean loggingOn) {

    instance->fmiVersion = FMIVersion2;

    instance->fmi2Functions = calloc(1, sizeof(FMI2Functions));

    if (!instance->fmi2Functions) {
        return FMIError;
    }

    instance->fmi2Functions->eventInfo.newDiscreteStatesNeeded           = fmi2False;
    instance->fmi2Functions->eventInfo.terminateSimulation               = fmi2False;
    instance->fmi2Functions->eventInfo.nominalsOfContinuousStatesChanged = fmi2False;
    instance->fmi2Functions->eventInfo.valuesOfContinuousStatesChanged   = fmi2False;
    instance->fmi2Functions->eventInfo.nextEventTimeDefined              = fmi2False;
    instance->fmi2Functions->eventInfo.nextEventTime                     = 0.0;

    instance->state = FMI2StartAndEndState;

#if !defined(FMI_VERSION) || FMI_VERSION == 2

    /***************************************************
    Common Functions
    ****************************************************/

    /* required functions */
    LOAD_SYMBOL(GetTypesPlatform);
    LOAD_SYMBOL(GetVersion);
    LOAD_SYMBOL(SetDebugLogging);
    LOAD_SYMBOL(Instantiate);
    LOAD_SYMBOL(FreeInstance);
    LOAD_SYMBOL(SetupExperiment);
    LOAD_SYMBOL(EnterInitializationMode);
    LOAD_SYMBOL(ExitInitializationMode);
    LOAD_SYMBOL(Terminate);
    LOAD_SYMBOL(Reset);
    LOAD_SYMBOL(GetReal);
    LOAD_SYMBOL(GetInteger);
    LOAD_SYMBOL(GetBoolean);
    LOAD_SYMBOL(GetString);
    LOAD_SYMBOL(SetReal);
    LOAD_SYMBOL(SetInteger);
    LOAD_SYMBOL(SetBoolean);
    LOAD_SYMBOL(SetString);

    /* optional functions */
    LOAD_SYMBOL(GetFMUstate);
    LOAD_SYMBOL(SetFMUstate);
    LOAD_SYMBOL(FreeFMUstate);
    LOAD_SYMBOL(SerializedFMUstateSize);
    LOAD_SYMBOL(SerializeFMUstate);
    LOAD_SYMBOL(DeSerializeFMUstate);
    LOAD_SYMBOL(GetDirectionalDerivative);

    if (fmuType == fmi2ModelExchange) {
#ifndef CO_SIMULATION
        /***************************************************
        Model Exchange
        ****************************************************/

        LOAD_SYMBOL(EnterEventMode);
        LOAD_SYMBOL(NewDiscreteStates);
        LOAD_SYMBOL(EnterContinuousTimeMode);
        LOAD_SYMBOL(CompletedIntegratorStep);
        LOAD_SYMBOL(SetTime);
        LOAD_SYMBOL(SetContinuousStates);
        LOAD_SYMBOL(GetDerivatives);
        LOAD_SYMBOL(GetEventIndicators);
        LOAD_SYMBOL(GetContinuousStates);
        LOAD_SYMBOL(GetNominalsOfContinuousStates);
#endif
    } else {
#ifndef MODEL_EXCHANGE
        /***************************************************
        Co-Simulation
        ****************************************************/

        LOAD_SYMBOL(SetRealInputDerivatives);
        LOAD_SYMBOL(GetRealOutputDerivatives);
        LOAD_SYMBOL(DoStep);
        LOAD_SYMBOL(CancelStep);
        LOAD_SYMBOL(GetStatus);
        LOAD_SYMBOL(GetRealStatus);
        LOAD_SYMBOL(GetIntegerStatus);
        LOAD_SYMBOL(GetBooleanStatus);
        LOAD_SYMBOL(GetStringStatus);
#endif
    }

#endif

    instance->fmi2Functions->callbacks.logger               = cb_logMessage2;
    instance->fmi2Functions->callbacks.allocateMemory       = calloc;
    instance->fmi2Functions->callbacks.freeMemory           = free;
    instance->fmi2Functions->callbacks.stepFinished         = NULL;
    instance->fmi2Functions->callbacks.componentEnvironment = instance;

    instance->component = instance->fmi2Functions->fmi2Instantiate(instance->name, fmuType, fmuGUID, fmuResourceLocation, &instance->fmi2Functions->callbacks, visible, loggingOn);

    if (instance->logFunctionCall) {
        fmi2CallbackFunctions* f = &instance->fmi2Functions->callbacks;
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi2Instantiate(instanceName=\"%s\", fmuType=%d, fmuGUID=\"%s\", fmuResourceLocation=\"%s\", functions={logger=0x%p, allocateMemory=0x%p, freeMemory=0x%p, stepFinished=0x%p, componentEnvironment=0x%p}, visible=%d, loggingOn=%d)",
            instance->name, fmuType, fmuGUID, fmuResourceLocation, f->logger, f->allocateMemory, f->freeMemory, f->stepFinished, f->componentEnvironment, visible, loggingOn);
        instance->logFunctionCall(instance, instance->component ? FMIOK : FMIError, instance->logMessageBuffer);
    }

    if (!instance->component) {
        return FMIError;
    }

    instance->interfaceType = (FMIInterfaceType)fmuType;
    instance->state = FMI2InstantiatedState;

    return FMIOK;
}

void FMI2FreeInstance(FMIInstance *instance) {

    if (!instance) {
        return;
    }

    instance->fmi2Functions->fmi2FreeInstance(instance->component);

    if (instance->logFunctionCall) {
        instance->logFunctionCall(instance, FMIOK, "fmi2FreeInstance()");
    }
}

/* Enter and exit initialization mode, terminate and reset */
FMIStatus FMI2SetupExperiment(FMIInstance *instance,
    fmi2Boolean toleranceDefined,
    fmi2Real tolerance,
    fmi2Real startTime,
    fmi2Boolean stopTimeDefined,
    fmi2Real stopTime) {

    instance->time = startTime;

    CALL_ARGS(SetupExperiment, "toleranceDefined=%d, tolerance=%.16g, startTime=%.16g, stopTimeDefined=%d, stopTime=%.16g", toleranceDefined, tolerance, startTime, stopTimeDefined, stopTime);
}

FMIStatus FMI2EnterInitializationMode(FMIInstance *instance) {
    instance->state = FMI2InitializationModeState;
    CALL(EnterInitializationMode);
}

FMIStatus FMI2ExitInitializationMode(FMIInstance *instance) {
    instance->state = instance->interfaceType == FMIModelExchange ? FMI2EventModeState : FMI2StepCompleteState;
    CALL(ExitInitializationMode);
}

FMIStatus FMI2Terminate(FMIInstance *instance) {
    instance->state = FMI2TerminatedState;
    CALL(Terminate);
}

FMIStatus FMI2Reset(FMIInstance *instance) {
    instance->state = FMI2InstantiatedState;
    CALL(Reset);
}

/* Getting and setting variable values */
FMIStatus FMI2GetReal(FMIInstance *instance, const fmi2ValueReference vr[], size_t nvr, fmi2Real value[]) {
    CALL_ARRAY(Get, Real);
}

FMIStatus FMI2GetInteger(FMIInstance *instance, const fmi2ValueReference vr[], size_t nvr, fmi2Integer value[]) {
    CALL_ARRAY(Get, Integer);
}

FMIStatus FMI2GetBoolean(FMIInstance *instance, const fmi2ValueReference vr[], size_t nvr, fmi2Boolean value[]) {
    CALL_ARRAY(Get, Boolean);
}

FMIStatus FMI2GetString(FMIInstance *instance, const fmi2ValueReference vr[], size_t nvr, fmi2String value[]) {
    CALL_ARRAY(Get, String);
}

FMIStatus FMI2SetReal(FMIInstance *instance, const fmi2ValueReference vr[], size_t nvr, const fmi2Real value[]) {
    CALL_ARRAY(Set, Real);
}

FMIStatus FMI2SetInteger(FMIInstance *instance, const fmi2ValueReference vr[], size_t nvr, const fmi2Integer value[]) {
    CALL_ARRAY(Set, Integer);
}

FMIStatus FMI2SetBoolean(FMIInstance *instance, const fmi2ValueReference vr[], size_t nvr, const fmi2Boolean value[]) {
    CALL_ARRAY(Set, Boolean);
}

FMIStatus FMI2SetString(FMIInstance *instance, const fmi2ValueReference vr[], size_t nvr, const fmi2String value[]) {
    CALL_ARRAY(Set, String);
}

/* Getting and setting the internal FMU state */
FMIStatus FMI2GetFMUstate(FMIInstance *instance, fmi2FMUstate* FMUstate) {
    CALL_ARGS(GetFMUstate, "FMUstate=0x%p", FMUstate);
}

FMIStatus FMI2SetFMUstate(FMIInstance *instance, fmi2FMUstate  FMUstate) {
    CALL_ARGS(SetFMUstate, "FMUstate=0x%p", FMUstate);
}

FMIStatus FMI2FreeFMUstate(FMIInstance *instance, fmi2FMUstate* FMUstate) {
    CALL_ARGS(FreeFMUstate, "FMUstate=0x%p", FMUstate);
}

FMIStatus FMI2SerializedFMUstateSize(FMIInstance *instance, fmi2FMUstate  FMUstate, size_t* size) {

    const FMIStatus status = (FMIStatus)instance->fmi2Functions->fmi2SerializedFMUstateSize(instance->component, FMUstate, size);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi2SerializedFMUstateSize(FMUstate=0x%p, size=%zu)", FMUstate, *size);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI2SerializeFMUstate(FMIInstance *instance, fmi2FMUstate  FMUstate, fmi2Byte serializedState[], size_t size) {
    CALL_ARGS(SerializeFMUstate, "FMUstate=0x%p, serializedState=0x%p, size=%zu", FMUstate, serializedState, size);
}

FMIStatus FMI2DeSerializeFMUstate(FMIInstance *instance, const fmi2Byte serializedState[], size_t size, fmi2FMUstate* FMUstate) {
    CALL_ARGS(DeSerializeFMUstate, "serializedState=0x%p, size=%zu, FMUstate=0x%p", serializedState, size, FMUstate);
}

/* Getting partial derivatives */
FMIStatus FMI2GetDirectionalDerivative(FMIInstance *instance,
    const fmi2ValueReference vUnknown_ref[], size_t nUnknown,
    const fmi2ValueReference vKnown_ref[], size_t nKnown,
    const fmi2Real dvKnown[],
    fmi2Real dvUnknown[]) {

    const FMIStatus status = (FMIStatus)instance->fmi2Functions->fmi2GetDirectionalDerivative(instance->component, vUnknown_ref, nUnknown, vKnown_ref, nKnown, dvKnown, dvUnknown);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi2GetDirectionalDerivative(vUnknown_ref={");
        FMIAppendArrayToLogMessageBuffer(instance, vUnknown_ref, nUnknown, NULL, FMIValueReferenceType);
        FMIAppendToLogMessageBuffer(instance, "}, nUnknown=%zu, vKnown_ref={", nUnknown);
        FMIAppendArrayToLogMessageBuffer(instance, vKnown_ref, nKnown, NULL, FMIValueReferenceType);
        FMIAppendToLogMessageBuffer(instance, "}, nKnown=%zu, dvKnown={", nKnown);
        FMIAppendArrayToLogMessageBuffer(instance, dvKnown, nKnown, NULL, FMIRealType);
        FMIAppendToLogMessageBuffer(instance, "}, dvUnknown={", nKnown);
        FMIAppendArrayToLogMessageBuffer(instance, dvUnknown, nUnknown, NULL, FMIRealType);
        FMIAppendToLogMessageBuffer(instance, "})");
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

/***************************************************
Model Exchange
****************************************************/

/* Enter and exit the different modes */
FMIStatus FMI2EnterEventMode(FMIInstance *instance) {
    instance->state = FMI2EventModeState;
    CALL(EnterEventMode);
}

FMIStatus FMI2NewDiscreteStates(FMIInstance *instance, fmi2EventInfo *eventInfo) {

    const FMIStatus status = (FMIStatus)instance->fmi2Functions->fmi2NewDiscreteStates(instance->component, eventInfo);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance,
            "fmi2NewDiscreteStates(eventInfo={newDiscreteStatesNeeded=%d, terminateSimulation=%d, nominalsOfContinuousStatesChanged=%d, valuesOfContinuousStatesChanged=%d, nextEventTimeDefined=%d, nextEventTime=%.16g})",
            eventInfo->newDiscreteStatesNeeded, eventInfo->terminateSimulation, eventInfo->nominalsOfContinuousStatesChanged, eventInfo->valuesOfContinuousStatesChanged, eventInfo->nextEventTimeDefined, eventInfo->nextEventTime);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI2EnterContinuousTimeMode(FMIInstance *instance) {
    instance->state = FMI2ContinuousTimeModeState;
    CALL(EnterContinuousTimeMode);
}

FMIStatus FMI2CompletedIntegratorStep(FMIInstance *instance,
    fmi2Boolean   noSetFMUStatePriorToCurrentPoint,
    fmi2Boolean*  enterEventMode,
    fmi2Boolean*  terminateSimulation) {

    const FMIStatus status = (FMIStatus)instance->fmi2Functions->fmi2CompletedIntegratorStep(instance->component, noSetFMUStatePriorToCurrentPoint, enterEventMode, terminateSimulation);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance,
            "fmi2CompletedIntegratorStep(noSetFMUStatePriorToCurrentPoint=%d, enterEventMode=%d, terminateSimulation=%d)",
            noSetFMUStatePriorToCurrentPoint, *enterEventMode, *terminateSimulation);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

/* Providing independent variables and re-initialization of caching */
FMIStatus FMI2SetTime(FMIInstance* instance, fmi2Real time) {
    instance->time = time;
    CALL_ARGS(SetTime, "time=%.16g", time);
}

FMIStatus FMI2SetContinuousStates(FMIInstance *instance, const fmi2Real x[], size_t nx) {

    const FMIStatus status = (FMIStatus)instance->fmi2Functions->fmi2SetContinuousStates(instance->component, x, nx);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi2SetContinuousStates(x={");
        FMIAppendArrayToLogMessageBuffer(instance, x, nx, NULL, FMIRealType);
        FMIAppendToLogMessageBuffer(instance, "}, nx=%zu)", nx);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

/* Evaluation of the model equations */
FMIStatus FMI2GetDerivatives(FMIInstance *instance, fmi2Real derivatives[], size_t nx) {

    const FMIStatus status = (FMIStatus)instance->fmi2Functions->fmi2GetDerivatives(instance->component, derivatives, nx);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi2GetDerivatives(derivatives={");
        FMIAppendArrayToLogMessageBuffer(instance, derivatives, nx, NULL, FMIRealType);
        FMIAppendToLogMessageBuffer(instance, "}, nx=%zu)", nx);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI2GetEventIndicators(FMIInstance *instance, fmi2Real eventIndicators[], size_t ni) {

    const FMIStatus status = (FMIStatus)instance->fmi2Functions->fmi2GetEventIndicators(instance->component, eventIndicators, ni);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi2GetEventIndicators(eventIndicators={");
        FMIAppendArrayToLogMessageBuffer(instance, eventIndicators, ni, NULL, FMIRealType);
        FMIAppendToLogMessageBuffer(instance, "}, ni=%zu)", ni);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI2GetContinuousStates(FMIInstance *instance, fmi2Real x[], size_t nx) {

    const FMIStatus status = (FMIStatus)instance->fmi2Functions->fmi2GetContinuousStates(instance->component, x, nx);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi2GetContinuousStates(x={");
        FMIAppendArrayToLogMessageBuffer(instance, x, nx, NULL, FMIRealType);
        FMIAppendToLogMessageBuffer(instance, "}, nx=%zu)", nx);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI2GetNominalsOfContinuousStates(FMIInstance *instance, fmi2Real x_nominal[], size_t nx) {

    const FMIStatus status = (FMIStatus)instance->fmi2Functions->fmi2GetNominalsOfContinuousStates(instance->component, x_nominal, nx);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi2GetNominalsOfContinuousStates(x_nominal={");
        FMIAppendArrayToLogMessageBuffer(instance, x_nominal, nx, NULL, FMIRealType);
        FMIAppendToLogMessageBuffer(instance, "}, nx=%zu)", nx);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

/***************************************************
Co-Simulation
****************************************************/

/* Simulating the slave */
FMIStatus FMI2SetRealInputDerivatives(FMIInstance *instance,
    const fmi2ValueReference vr[], size_t nvr,
    const fmi2Integer order[],
    const fmi2Real value[]) {
    CALL_ARGS(SetRealInputDerivatives, "vr=0x%p, nvr=%zu, order=0x%p, value=0x%p", vr, nvr, order, value);
}

FMIStatus FMI2GetRealOutputDerivatives(FMIInstance *instance,
    const fmi2ValueReference vr[], size_t nvr,
    const fmi2Integer order[],
    fmi2Real value[]) {
    CALL_ARGS(GetRealOutputDerivatives, "vr=0x%p, nvr=%zu, order=0x%p, value=0x%p", vr, nvr, order, value);
}

FMIStatus FMI2DoStep(FMIInstance *instance,
    fmi2Real      currentCommunicationPoint,
    fmi2Real      communicationStepSize,
    fmi2Boolean   noSetFMUStatePriorToCurrentPoint) {

    instance->time = currentCommunicationPoint + communicationStepSize;

    CALL_ARGS(DoStep, "currentCommunicationPoint=%.16g, communicationStepSize=%.16g, noSetFMUStatePriorToCurrentPoint=%d",
        currentCommunicationPoint, communicationStepSize, noSetFMUStatePriorToCurrentPoint);
}

FMIStatus FMI2CancelStep(FMIInstance *instance) {
    CALL(CancelStep);
}

/* Inquire slave status */
FMIStatus FMI2GetStatus(FMIInstance *instance, const fmi2StatusKind s, fmi2Status* value) {

    const FMIStatus status = (FMIStatus)instance->fmi2Functions->fmi2GetStatus(instance->component, s, value);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi2GetStatus(s=%d, value=%d)", s, *value);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI2GetRealStatus(FMIInstance *instance, const fmi2StatusKind s, fmi2Real* value) {

    const FMIStatus status = (FMIStatus)instance->fmi2Functions->fmi2GetRealStatus(instance->component, s, value);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi2GetRealStatus(s=%d, value=%.16g)", s, *value);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI2GetIntegerStatus(FMIInstance *instance, const fmi2StatusKind s, fmi2Integer* value) {

    const FMIStatus status = (FMIStatus)instance->fmi2Functions->fmi2GetIntegerStatus(instance->component, s, value);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi2GetIntegerStatus(s=%d, value=%d)", s, *value);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI2GetBooleanStatus(FMIInstance *instance, const fmi2StatusKind s, fmi2Boolean* value) {

    const FMIStatus status = (FMIStatus)instance->fmi2Functions->fmi2GetBooleanStatus(instance->component, s, value);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi2GetBooleanStatus(s=%d, value=%d)", s, *value);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI2GetStringStatus(FMIInstance *instance, const fmi2StatusKind s, fmi2String* value) {

    const FMIStatus status = (FMIStatus)instance->fmi2Functions->fmi2GetStringStatus(instance->component, s, value);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi2GetStringStatus(s=%d, value=\"%s\")", s, *value);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

#undef LOAD_SYMBOL
#undef CALL
#undef CALL_ARGS
#undef CALL_ARRAY
#ifdef _WIN32
#include <direct.h>
#include "Shlwapi.h"
#pragma comment(lib, "shlwapi.lib")
#else
#include <stdarg.h>
#include <dlfcn.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#ifndef fmi3FunctionTypes_h
#define fmi3FunctionTypes_h

#ifndef fmi3PlatformTypes_h
#define fmi3PlatformTypes_h

/*
This header file defines the data types of FMI 3.0.
It must be used by both FMU and FMI master.

Copyright (C) 2008-2011 MODELISAR consortium,
              2012-2022 Modelica Association Project "FMI"
              All rights reserved.

This file is licensed by the copyright holders under the 2-Clause BSD License
(https://opensource.org/licenses/BSD-2-Clause):

----------------------------------------------------------------------------
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
----------------------------------------------------------------------------
*/

/* Include the integer and boolean type definitions */
#include <stdint.h>
#include <stdbool.h>


/* tag::Component[] */
typedef           void* fmi3Instance;             /* Pointer to the FMU instance */
/* end::Component[] */

/* tag::ComponentEnvironment[] */
typedef           void* fmi3InstanceEnvironment;  /* Pointer to the FMU environment */
/* end::ComponentEnvironment[] */

/* tag::FMUState[] */
typedef           void* fmi3FMUState;             /* Pointer to the internal FMU state */
/* end::FMUState[] */

/* tag::ValueReference[] */
typedef        uint32_t fmi3ValueReference;       /* Handle to the value of a variable */
/* end::ValueReference[] */

/* tag::VariableTypes[] */
typedef           float fmi3Float32;  /* Single precision floating point (32-bit) */
/* tag::fmi3Float64[] */
typedef          double fmi3Float64;  /* Double precision floating point (64-bit) */
/* end::fmi3Float64[] */
typedef          int8_t fmi3Int8;     /* 8-bit signed integer */
typedef         uint8_t fmi3UInt8;    /* 8-bit unsigned integer */
typedef         int16_t fmi3Int16;    /* 16-bit signed integer */
typedef        uint16_t fmi3UInt16;   /* 16-bit unsigned integer */
typedef         int32_t fmi3Int32;    /* 32-bit signed integer */
typedef        uint32_t fmi3UInt32;   /* 32-bit unsigned integer */
typedef         int64_t fmi3Int64;    /* 64-bit signed integer */
typedef        uint64_t fmi3UInt64;   /* 64-bit unsigned integer */
typedef            bool fmi3Boolean;  /* Data type to be used with fmi3True and fmi3False */
typedef            char fmi3Char;     /* Data type for one character */
typedef const fmi3Char* fmi3String;   /* Data type for character strings
                                         ('\0' terminated, UTF-8 encoded) */
typedef         uint8_t fmi3Byte;     /* Smallest addressable unit of the machine
                                         (typically one byte) */
typedef const fmi3Byte* fmi3Binary;   /* Data type for binary data
                                         (out-of-band length terminated) */
typedef            bool fmi3Clock;    /* Data type to be used with fmi3ClockActive and
                                         fmi3ClockInactive */

/* Values for fmi3Boolean */
#define fmi3True  true
#define fmi3False false

/* Values for fmi3Clock */
#define fmi3ClockActive   true
#define fmi3ClockInactive false
/* end::VariableTypes[] */

#endif /* fmi3PlatformTypes_h */

/*
This header file defines the data and function types of FMI 3.0.
It must be used when compiling an FMU or an FMI importer.

Copyright (C) 2011 MODELISAR consortium,
              2012-2022 Modelica Association Project "FMI"
              All rights reserved.

This file is licensed by the copyright holders under the 2-Clause BSD License
(https://opensource.org/licenses/BSD-2-Clause):

----------------------------------------------------------------------------
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
----------------------------------------------------------------------------
*/

#ifdef __cplusplus
extern "C" {
#endif

/* Include stddef.h, in order that size_t etc. is defined */
#include <stddef.h>


/* Type definitions */

/* tag::Status[] */
typedef enum {
    fmi3OK,
    fmi3Warning,
    fmi3Discard,
    fmi3Error,
    fmi3Fatal,
} fmi3Status;
/* end::Status[] */

/* tag::DependencyKind[] */
typedef enum {
    fmi3Independent,
    fmi3Constant,
    fmi3Fixed,
    fmi3Tunable,
    fmi3Discrete,
    fmi3Dependent
} fmi3DependencyKind;
/* end::DependencyKind[] */

/* tag::IntervalQualifier[] */
typedef enum {
    fmi3IntervalNotYetKnown,
    fmi3IntervalUnchanged,
    fmi3IntervalChanged
} fmi3IntervalQualifier;
/* end::IntervalQualifier[] */

/* tag::CallbackLogMessage[] */
typedef void  (*fmi3LogMessageCallback) (fmi3InstanceEnvironment instanceEnvironment,
                                         fmi3Status status,
                                         fmi3String category,
                                         fmi3String message);
/* end::CallbackLogMessage[] */

/* tag::CallbackClockUpdate[] */
typedef void (*fmi3ClockUpdateCallback) (
    fmi3InstanceEnvironment  instanceEnvironment);
/* end::CallbackClockUpdate[] */

/* tag::CallbackIntermediateUpdate[] */
typedef void (*fmi3IntermediateUpdateCallback) (
    fmi3InstanceEnvironment instanceEnvironment,
    fmi3Float64  intermediateUpdateTime,
    fmi3Boolean  intermediateVariableSetRequested,
    fmi3Boolean  intermediateVariableGetAllowed,
    fmi3Boolean  intermediateStepFinished,
    fmi3Boolean  canReturnEarly,
    fmi3Boolean* earlyReturnRequested,
    fmi3Float64* earlyReturnTime);
/* end::CallbackIntermediateUpdate[] */

/* tag::CallbackPreemptionLock[] */
typedef void (*fmi3LockPreemptionCallback)   ();
typedef void (*fmi3UnlockPreemptionCallback) ();
/* end::CallbackPreemptionLock[] */

/* Define fmi3 function pointer types to simplify dynamic loading */

/***************************************************
Types for Common Functions
****************************************************/

/* Inquire version numbers and setting logging status */
/* tag::GetVersion[] */
typedef const char* fmi3GetVersionTYPE(void);
/* end::GetVersion[] */

/* tag::SetDebugLogging[] */
typedef fmi3Status fmi3SetDebugLoggingTYPE(fmi3Instance instance,
                                           fmi3Boolean loggingOn,
                                           size_t nCategories,
                                           const fmi3String categories[]);
/* end::SetDebugLogging[] */

/* Creation and destruction of FMU instances and setting debug status */
/* tag::Instantiate[] */
typedef fmi3Instance fmi3InstantiateModelExchangeTYPE(
    fmi3String                 instanceName,
    fmi3String                 instantiationToken,
    fmi3String                 resourcePath,
    fmi3Boolean                visible,
    fmi3Boolean                loggingOn,
    fmi3InstanceEnvironment    instanceEnvironment,
    fmi3LogMessageCallback     logMessage);

typedef fmi3Instance fmi3InstantiateCoSimulationTYPE(
    fmi3String                     instanceName,
    fmi3String                     instantiationToken,
    fmi3String                     resourcePath,
    fmi3Boolean                    visible,
    fmi3Boolean                    loggingOn,
    fmi3Boolean                    eventModeUsed,
    fmi3Boolean                    earlyReturnAllowed,
    const fmi3ValueReference       requiredIntermediateVariables[],
    size_t                         nRequiredIntermediateVariables,
    fmi3InstanceEnvironment        instanceEnvironment,
    fmi3LogMessageCallback         logMessage,
    fmi3IntermediateUpdateCallback intermediateUpdate);

typedef fmi3Instance fmi3InstantiateScheduledExecutionTYPE(
    fmi3String                     instanceName,
    fmi3String                     instantiationToken,
    fmi3String                     resourcePath,
    fmi3Boolean                    visible,
    fmi3Boolean                    loggingOn,
    fmi3InstanceEnvironment        instanceEnvironment,
    fmi3LogMessageCallback         logMessage,
    fmi3ClockUpdateCallback        clockUpdate,
    fmi3LockPreemptionCallback     lockPreemption,
    fmi3UnlockPreemptionCallback   unlockPreemption);
/* end::Instantiate[] */

/* tag::FreeInstance[] */
typedef void fmi3FreeInstanceTYPE(fmi3Instance instance);
/* end::FreeInstance[] */

/* Enter and exit initialization mode, enter event mode, terminate and reset */
/* tag::EnterInitializationMode[] */
typedef fmi3Status fmi3EnterInitializationModeTYPE(fmi3Instance instance,
                                                   fmi3Boolean toleranceDefined,
                                                   fmi3Float64 tolerance,
                                                   fmi3Float64 startTime,
                                                   fmi3Boolean stopTimeDefined,
                                                   fmi3Float64 stopTime);
/* end::EnterInitializationMode[] */

/* tag::ExitInitializationMode[] */
typedef fmi3Status fmi3ExitInitializationModeTYPE(fmi3Instance instance);
/* end::ExitInitializationMode[] */

/* tag::EnterEventMode[] */
typedef fmi3Status fmi3EnterEventModeTYPE(fmi3Instance instance);
/* end::EnterEventMode[] */

/* tag::Terminate[] */
typedef fmi3Status fmi3TerminateTYPE(fmi3Instance instance);
/* end::Terminate[] */

/* tag::Reset[] */
typedef fmi3Status fmi3ResetTYPE(fmi3Instance instance);
/* end::Reset[] */

/* Getting and setting variable values */
/* tag::Getters[] */
typedef fmi3Status fmi3GetFloat32TYPE(fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3Float32 values[],
                                      size_t nValues);

typedef fmi3Status fmi3GetFloat64TYPE(fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3Float64 values[],
                                      size_t nValues);

typedef fmi3Status fmi3GetInt8TYPE   (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3Int8 values[],
                                      size_t nValues);

typedef fmi3Status fmi3GetUInt8TYPE  (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3UInt8 values[],
                                      size_t nValues);

typedef fmi3Status fmi3GetInt16TYPE  (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3Int16 values[],
                                      size_t nValues);

typedef fmi3Status fmi3GetUInt16TYPE (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3UInt16 values[],
                                      size_t nValues);

typedef fmi3Status fmi3GetInt32TYPE  (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3Int32 values[],
                                      size_t nValues);

typedef fmi3Status fmi3GetUInt32TYPE (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3UInt32 values[],
                                      size_t nValues);

typedef fmi3Status fmi3GetInt64TYPE  (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3Int64 values[],
                                      size_t nValues);

typedef fmi3Status fmi3GetUInt64TYPE (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3UInt64 values[],
                                      size_t nValues);

typedef fmi3Status fmi3GetBooleanTYPE(fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3Boolean values[],
                                      size_t nValues);

typedef fmi3Status fmi3GetStringTYPE (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3String values[],
                                      size_t nValues);

typedef fmi3Status fmi3GetBinaryTYPE (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      size_t valueSizes[],
                                      fmi3Binary values[],
                                      size_t nValues);
/* end::Getters[] */

/* tag::GetClock[] */
typedef fmi3Status fmi3GetClockTYPE  (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3Clock values[]);
/* end::GetClock[] */

/* tag::Setters[] */
typedef fmi3Status fmi3SetFloat32TYPE(fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3Float32 values[],
                                      size_t nValues);

typedef fmi3Status fmi3SetFloat64TYPE(fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3Float64 values[],
                                      size_t nValues);

typedef fmi3Status fmi3SetInt8TYPE   (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3Int8 values[],
                                      size_t nValues);

typedef fmi3Status fmi3SetUInt8TYPE  (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3UInt8 values[],
                                      size_t nValues);

typedef fmi3Status fmi3SetInt16TYPE  (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3Int16 values[],
                                      size_t nValues);

typedef fmi3Status fmi3SetUInt16TYPE (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3UInt16 values[],
                                      size_t nValues);

typedef fmi3Status fmi3SetInt32TYPE  (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3Int32 values[],
                                      size_t nValues);

typedef fmi3Status fmi3SetUInt32TYPE (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3UInt32 values[],
                                      size_t nValues);

typedef fmi3Status fmi3SetInt64TYPE  (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3Int64 values[],
                                      size_t nValues);

typedef fmi3Status fmi3SetUInt64TYPE (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3UInt64 values[],
                                      size_t nValues);

typedef fmi3Status fmi3SetBooleanTYPE(fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3Boolean values[],
                                      size_t nValues);

typedef fmi3Status fmi3SetStringTYPE (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3String values[],
                                      size_t nValues);

typedef fmi3Status fmi3SetBinaryTYPE (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const size_t valueSizes[],
                                      const fmi3Binary values[],
                                      size_t nValues);
/* end::Setters[] */
/* tag::SetClock[] */
typedef fmi3Status fmi3SetClockTYPE  (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3Clock values[]);
/* end::SetClock[] */

/* Getting Variable Dependency Information */
/* tag::GetNumberOfVariableDependencies[] */
typedef fmi3Status fmi3GetNumberOfVariableDependenciesTYPE(fmi3Instance instance,
                                                           fmi3ValueReference valueReference,
                                                           size_t* nDependencies);
/* end::GetNumberOfVariableDependencies[] */

/* tag::GetVariableDependencies[] */
typedef fmi3Status fmi3GetVariableDependenciesTYPE(fmi3Instance instance,
                                                   fmi3ValueReference dependent,
                                                   size_t elementIndicesOfDependent[],
                                                   fmi3ValueReference independents[],
                                                   size_t elementIndicesOfIndependents[],
                                                   fmi3DependencyKind dependencyKinds[],
                                                   size_t nDependencies);
/* end::GetVariableDependencies[] */

/* Getting and setting the internal FMU state */
/* tag::GetFMUState[] */
typedef fmi3Status fmi3GetFMUStateTYPE (fmi3Instance instance, fmi3FMUState* FMUState);
/* end::GetFMUState[] */

/* tag::SetFMUState[] */
typedef fmi3Status fmi3SetFMUStateTYPE (fmi3Instance instance, fmi3FMUState  FMUState);
/* end::SetFMUState[] */

/* tag::FreeFMUState[] */
typedef fmi3Status fmi3FreeFMUStateTYPE(fmi3Instance instance, fmi3FMUState* FMUState);
/* end::FreeFMUState[] */

/* tag::SerializedFMUStateSize[] */
typedef fmi3Status fmi3SerializedFMUStateSizeTYPE(fmi3Instance instance,
                                                  fmi3FMUState FMUState,
                                                  size_t* size);
/* end::SerializedFMUStateSize[] */

/* tag::SerializeFMUState[] */
typedef fmi3Status fmi3SerializeFMUStateTYPE     (fmi3Instance instance,
                                                  fmi3FMUState FMUState,
                                                  fmi3Byte serializedState[],
                                                  size_t size);
/* end::SerializeFMUState[] */

/* tag::DeserializeFMUState[] */
typedef fmi3Status fmi3DeserializeFMUStateTYPE   (fmi3Instance instance,
                                                  const fmi3Byte serializedState[],
                                                  size_t size,
                                                  fmi3FMUState* FMUState);
/* end::DeserializeFMUState[] */

/* Getting partial derivatives */
/* tag::GetDirectionalDerivative[] */
typedef fmi3Status fmi3GetDirectionalDerivativeTYPE(fmi3Instance instance,
                                                    const fmi3ValueReference unknowns[],
                                                    size_t nUnknowns,
                                                    const fmi3ValueReference knowns[],
                                                    size_t nKnowns,
                                                    const fmi3Float64 seed[],
                                                    size_t nSeed,
                                                    fmi3Float64 sensitivity[],
                                                    size_t nSensitivity);
/* end::GetDirectionalDerivative[] */

/* tag::GetAdjointDerivative[] */
typedef fmi3Status fmi3GetAdjointDerivativeTYPE(fmi3Instance instance,
                                                const fmi3ValueReference unknowns[],
                                                size_t nUnknowns,
                                                const fmi3ValueReference knowns[],
                                                size_t nKnowns,
                                                const fmi3Float64 seed[],
                                                size_t nSeed,
                                                fmi3Float64 sensitivity[],
                                                size_t nSensitivity);
/* end::GetAdjointDerivative[] */

/* Entering and exiting the Configuration or Reconfiguration Mode */

/* tag::EnterConfigurationMode[] */
typedef fmi3Status fmi3EnterConfigurationModeTYPE(fmi3Instance instance);
/* end::EnterConfigurationMode[] */

/* tag::ExitConfigurationMode[] */
typedef fmi3Status fmi3ExitConfigurationModeTYPE(fmi3Instance instance);
/* end::ExitConfigurationMode[] */

/* tag::GetIntervalDecimal[] */
typedef fmi3Status fmi3GetIntervalDecimalTYPE(fmi3Instance instance,
                                              const fmi3ValueReference valueReferences[],
                                              size_t nValueReferences,
                                              fmi3Float64 intervals[],
                                              fmi3IntervalQualifier qualifiers[]);
/* end::GetIntervalDecimal[] */

/* tag::GetIntervalFraction[] */
typedef fmi3Status fmi3GetIntervalFractionTYPE(fmi3Instance instance,
                                               const fmi3ValueReference valueReferences[],
                                               size_t nValueReferences,
                                               fmi3UInt64 counters[],
                                               fmi3UInt64 resolutions[],
                                               fmi3IntervalQualifier qualifiers[]);
/* end::GetIntervalFraction[] */

/* tag::GetShiftDecimal[] */
typedef fmi3Status fmi3GetShiftDecimalTYPE(fmi3Instance instance,
                                           const fmi3ValueReference valueReferences[],
                                           size_t nValueReferences,
                                           fmi3Float64 shifts[]);
/* end::GetShiftDecimal[] */

/* tag::GetShiftFraction[] */
typedef fmi3Status fmi3GetShiftFractionTYPE(fmi3Instance instance,
                                            const fmi3ValueReference valueReferences[],
                                            size_t nValueReferences,
                                            fmi3UInt64 counters[],
                                            fmi3UInt64 resolutions[]);
/* end::GetShiftFraction[] */

/* tag::SetIntervalDecimal[] */
typedef fmi3Status fmi3SetIntervalDecimalTYPE(fmi3Instance instance,
                                              const fmi3ValueReference valueReferences[],
                                              size_t nValueReferences,
                                              const fmi3Float64 intervals[]);
/* end::SetIntervalDecimal[] */

/* tag::SetIntervalFraction[] */
typedef fmi3Status fmi3SetIntervalFractionTYPE(fmi3Instance instance,
                                               const fmi3ValueReference valueReferences[],
                                               size_t nValueReferences,
                                               const fmi3UInt64 counters[],
                                               const fmi3UInt64 resolutions[]);
/* end::SetIntervalFraction[] */

/* tag::SetShiftDecimal[] */
typedef fmi3Status fmi3SetShiftDecimalTYPE(fmi3Instance instance,
                                           const fmi3ValueReference valueReferences[],
                                           size_t nValueReferences,
                                           const fmi3Float64 shifts[]);
/* end::SetShiftDecimal[] */

/* tag::SetShiftFraction[] */
typedef fmi3Status fmi3SetShiftFractionTYPE(fmi3Instance instance,
                                            const fmi3ValueReference valueReferences[],
                                            size_t nValueReferences,
                                            const fmi3UInt64 counters[],
                                            const fmi3UInt64 resolutions[]);
/* end::SetShiftFraction[] */

/* tag::EvaluateDiscreteStates[] */
typedef fmi3Status fmi3EvaluateDiscreteStatesTYPE(fmi3Instance instance);
/* end::EvaluateDiscreteStates[] */

/* tag::UpdateDiscreteStates[] */
typedef fmi3Status fmi3UpdateDiscreteStatesTYPE(fmi3Instance instance,
                                                fmi3Boolean* discreteStatesNeedUpdate,
                                                fmi3Boolean* terminateSimulation,
                                                fmi3Boolean* nominalsOfContinuousStatesChanged,
                                                fmi3Boolean* valuesOfContinuousStatesChanged,
                                                fmi3Boolean* nextEventTimeDefined,
                                                fmi3Float64* nextEventTime);
/* end::UpdateDiscreteStates[] */

/***************************************************
Types for Functions for Model Exchange
****************************************************/

/* tag::EnterContinuousTimeMode[] */
typedef fmi3Status fmi3EnterContinuousTimeModeTYPE(fmi3Instance instance);
/* end::EnterContinuousTimeMode[] */

/* tag::CompletedIntegratorStep[] */
typedef fmi3Status fmi3CompletedIntegratorStepTYPE(fmi3Instance instance,
                                                   fmi3Boolean  noSetFMUStatePriorToCurrentPoint,
                                                   fmi3Boolean* enterEventMode,
                                                   fmi3Boolean* terminateSimulation);
/* end::CompletedIntegratorStep[] */

/* Providing independent variables and re-initialization of caching */
/* tag::SetTime[] */
typedef fmi3Status fmi3SetTimeTYPE(fmi3Instance instance, fmi3Float64 time);
/* end::SetTime[] */

/* tag::SetContinuousStates[] */
typedef fmi3Status fmi3SetContinuousStatesTYPE(fmi3Instance instance,
                                               const fmi3Float64 continuousStates[],
                                               size_t nContinuousStates);
/* end::SetContinuousStates[] */

/* Evaluation of the model equations */
/* tag::GetDerivatives[] */
typedef fmi3Status fmi3GetContinuousStateDerivativesTYPE(fmi3Instance instance,
                                                         fmi3Float64 derivatives[],
                                                         size_t nContinuousStates);
/* end::GetDerivatives[] */

/* tag::GetEventIndicators[] */
typedef fmi3Status fmi3GetEventIndicatorsTYPE(fmi3Instance instance,
                                              fmi3Float64 eventIndicators[],
                                              size_t nEventIndicators);
/* end::GetEventIndicators[] */

/* tag::GetContinuousStates[] */
typedef fmi3Status fmi3GetContinuousStatesTYPE(fmi3Instance instance,
                                               fmi3Float64 continuousStates[],
                                               size_t nContinuousStates);
/* end::GetContinuousStates[] */

/* tag::GetNominalsOfContinuousStates[] */
typedef fmi3Status fmi3GetNominalsOfContinuousStatesTYPE(fmi3Instance instance,
                                                         fmi3Float64 nominals[],
                                                         size_t nContinuousStates);
/* end::GetNominalsOfContinuousStates[] */

/* tag::GetNumberOfEventIndicators[] */
typedef fmi3Status fmi3GetNumberOfEventIndicatorsTYPE(fmi3Instance instance,
                                                      size_t* nEventIndicators);
/* end::GetNumberOfEventIndicators[] */

/* tag::GetNumberOfContinuousStates[] */
typedef fmi3Status fmi3GetNumberOfContinuousStatesTYPE(fmi3Instance instance,
                                                       size_t* nContinuousStates);
/* end::GetNumberOfContinuousStates[] */

/***************************************************
Types for Functions for Co-Simulation
****************************************************/

/* Simulating the FMU */

/* tag::EnterStepMode[] */
typedef fmi3Status fmi3EnterStepModeTYPE(fmi3Instance instance);
/* end::EnterStepMode[] */

/* tag::GetOutputDerivatives[] */
typedef fmi3Status fmi3GetOutputDerivativesTYPE(fmi3Instance instance,
                                                const fmi3ValueReference valueReferences[],
                                                size_t nValueReferences,
                                                const fmi3Int32 orders[],
                                                fmi3Float64 values[],
                                                size_t nValues);
/* end::GetOutputDerivatives[] */

/* tag::DoStep[] */
typedef fmi3Status fmi3DoStepTYPE(fmi3Instance instance,
                                  fmi3Float64 currentCommunicationPoint,
                                  fmi3Float64 communicationStepSize,
                                  fmi3Boolean noSetFMUStatePriorToCurrentPoint,
                                  fmi3Boolean* eventHandlingNeeded,
                                  fmi3Boolean* terminateSimulation,
                                  fmi3Boolean* earlyReturn,
                                  fmi3Float64* lastSuccessfulTime);
/* end::DoStep[] */

/***************************************************
Types for Functions for Scheduled Execution
****************************************************/

/* tag::ActivateModelPartition[] */
typedef fmi3Status fmi3ActivateModelPartitionTYPE(fmi3Instance instance,
                                                  fmi3ValueReference clockReference,
                                                  fmi3Float64 activationTime);
/* end::ActivateModelPartition[] */

#ifdef __cplusplus
}  /* end of extern "C" { */
#endif

#endif /* fmi3FunctionTypes_h */
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <Windows.h>
#endif

#include <stdbool.h>
#include <stddef.h>

#ifndef FMI_MAX_MESSAGE_LENGTH
#define FMI_MAX_MESSAGE_LENGTH 4096
#endif

#ifndef FMI_STATIC
#define FMI_STATIC
#endif

typedef enum {
    FMIOK,
    FMIWarning,
    FMIDiscard,
    FMIError,
    FMIFatal,
    FMIPending
} FMIStatus;

typedef enum {

    // FMI 3.0 variable types
    FMIFloat32Type,
    FMIDiscreteFloat32Type,
    FMIFloat64Type,
    FMIDiscreteFloat64Type,
    FMIInt8Type,
    FMIUInt8Type,
    FMIInt16Type,
    FMIUInt16Type,
    FMIInt32Type,
    FMIUInt32Type,
    FMIInt64Type,
    FMIUInt64Type,
    FMIBooleanType,
    FMIStringType,
    FMIBinaryType,
    FMIClockType,

    // non-variable types used as parameters
    FMIValueReferenceType,
    FMISizeTType,

    // aliases for FMI 1.0 and 2.0 variable types
    FMIRealType         = FMIFloat64Type,
    FMIDiscreteRealType = FMIDiscreteFloat64Type,
    FMIIntegerType      = FMIInt32Type,

} FMIVariableType;

typedef enum {
    FMIVersion1 = 1,
    FMIVersion2 = 2,
    FMIVersion3 = 3
} FMIVersion;

typedef enum {
    FMIModelExchange,
    FMICoSimulation,
    FMIScheduledExecution
} FMIInterfaceType;

typedef enum {
    FMI2StartAndEndState        = 1 << 0,
    FMI2InstantiatedState       = 1 << 1,
    FMI2InitializationModeState = 1 << 2,

    // model exchange states
    FMI2EventModeState          = 1 << 3,
    FMI2ContinuousTimeModeState = 1 << 4,

    // co-simulation states
    FMI2StepCompleteState       = 1 << 5,
    FMI2StepInProgressState     = 1 << 6,
    FMI2StepFailedState         = 1 << 7,
    FMI2StepCanceledState       = 1 << 8,

    FMI2TerminatedState         = 1 << 9,
    FMI2ErrorState              = 1 << 10,
    FMI2FatalState              = 1 << 11,
} FMI2State;

typedef unsigned int FMIValueReference;

typedef struct FMIInstance_ FMIInstance;

typedef struct FMI1Functions_ FMI1Functions;

typedef struct FMI2Functions_ FMI2Functions;

typedef struct FMI3Functions_ FMI3Functions;

typedef void FMILogFunctionCall(FMIInstance *instance, FMIStatus status, const char *message);

typedef void FMILogMessage(FMIInstance *instance, FMIStatus status, const char *category, const char *message);

struct FMIInstance_ {

    FMI1Functions *fmi1Functions;
    FMI2Functions *fmi2Functions;
    FMI3Functions *fmi3Functions;

#ifdef _WIN32
    HMODULE libraryHandle;
#else
    void *libraryHandle;
#endif

    void *userData;

    FMILogMessage      *logMessage;
    FMILogFunctionCall *logFunctionCall;

    double time;

    char* logMessageBuffer;
    size_t logMessageBufferSize;
    size_t logMessageBufferPosition;

    void *component;

    const char *name;

    bool logFMICalls;

    FMI2State state;

    FMIStatus status;

    FMIVersion fmiVersion;

    FMIInterfaceType interfaceType;

};

FMI_STATIC FMIInstance* FMICreateInstance(const char* instanceName, const char* libraryPath, FMILogMessage* logMessage, FMILogFunctionCall* logFunctionCall);

FMI_STATIC void FMIFreeInstance(FMIInstance *instance);

FMI_STATIC void FMIClearLogMessageBuffer(FMIInstance* instance);

FMI_STATIC void FMIAppendToLogMessageBuffer(FMIInstance* instance, const char* format, ...);

FMI_STATIC void FMIAppendArrayToLogMessageBuffer(FMIInstance* instance, const void* values, size_t nValues, const size_t sizes[], FMIVariableType variableType);

FMI_STATIC FMIStatus FMIURIToPath(const char *uri, char *path, const size_t pathLength);

FMI_STATIC FMIStatus FMIPathToURI(const char *path, char *uri, const size_t uriLength);

FMI_STATIC FMIStatus FMIPlatformBinaryPath(const char *unzipdir, const char *modelIdentifier, FMIVersion fmiVersion, char *platformBinaryPath, size_t size);

#ifdef __cplusplus
}  /* end of extern "C" { */
#endif

struct FMI3Functions_ {

    fmi3Boolean eventModeUsed;

    fmi3Boolean discreteStatesNeedUpdate;
    fmi3Boolean terminateSimulation;
    fmi3Boolean nominalsOfContinuousStatesChanged;
    fmi3Boolean valuesOfContinuousStatesChanged;
    fmi3Boolean nextEventTimeDefined;
    fmi3Float64 nextEventTime;

    /***************************************************
    Common Functions for FMI 3.0
    ****************************************************/

    /* Inquire version numbers and set debug logging */
    fmi3GetVersionTYPE                      *fmi3GetVersion;
    fmi3SetDebugLoggingTYPE                 *fmi3SetDebugLogging;

    /* Creation and destruction of FMU instances */
    fmi3InstantiateModelExchangeTYPE        *fmi3InstantiateModelExchange;
    fmi3InstantiateCoSimulationTYPE         *fmi3InstantiateCoSimulation;
    fmi3InstantiateScheduledExecutionTYPE   *fmi3InstantiateScheduledExecution;
    fmi3FreeInstanceTYPE                    *fmi3FreeInstance;

    /* Enter and exit initialization mode, terminate and reset */
    fmi3EnterInitializationModeTYPE         *fmi3EnterInitializationMode;
    fmi3ExitInitializationModeTYPE          *fmi3ExitInitializationMode;
    fmi3EnterEventModeTYPE                  *fmi3EnterEventMode;
    fmi3TerminateTYPE                       *fmi3Terminate;
    fmi3ResetTYPE                           *fmi3Reset;

    /* Getting and setting variable values */
    fmi3GetFloat32TYPE                      *fmi3GetFloat32;
    fmi3GetFloat64TYPE                      *fmi3GetFloat64;
    fmi3GetInt8TYPE                         *fmi3GetInt8;
    fmi3GetUInt8TYPE                        *fmi3GetUInt8;
    fmi3GetInt16TYPE                        *fmi3GetInt16;
    fmi3GetUInt16TYPE                       *fmi3GetUInt16;
    fmi3GetInt32TYPE                        *fmi3GetInt32;
    fmi3GetUInt32TYPE                       *fmi3GetUInt32;
    fmi3GetInt64TYPE                        *fmi3GetInt64;
    fmi3GetUInt64TYPE                       *fmi3GetUInt64;
    fmi3GetBooleanTYPE                      *fmi3GetBoolean;
    fmi3GetStringTYPE                       *fmi3GetString;
    fmi3GetBinaryTYPE                       *fmi3GetBinary;
    fmi3GetClockTYPE                        *fmi3GetClock;
    fmi3SetFloat32TYPE                      *fmi3SetFloat32;
    fmi3SetFloat64TYPE                      *fmi3SetFloat64;
    fmi3SetInt8TYPE                         *fmi3SetInt8;
    fmi3SetUInt8TYPE                        *fmi3SetUInt8;
    fmi3SetInt16TYPE                        *fmi3SetInt16;
    fmi3SetUInt16TYPE                       *fmi3SetUInt16;
    fmi3SetInt32TYPE                        *fmi3SetInt32;
    fmi3SetUInt32TYPE                       *fmi3SetUInt32;
    fmi3SetInt64TYPE                        *fmi3SetInt64;
    fmi3SetUInt64TYPE                       *fmi3SetUInt64;
    fmi3SetBooleanTYPE                      *fmi3SetBoolean;
    fmi3SetStringTYPE                       *fmi3SetString;
    fmi3SetBinaryTYPE                       *fmi3SetBinary;
    fmi3SetClockTYPE                        *fmi3SetClock;

    /* Getting Variable Dependency Information */
    fmi3GetNumberOfVariableDependenciesTYPE *fmi3GetNumberOfVariableDependencies;
    fmi3GetVariableDependenciesTYPE         *fmi3GetVariableDependencies;

    /* Getting and setting the internal FMU state */
    fmi3GetFMUStateTYPE                     *fmi3GetFMUState;
    fmi3SetFMUStateTYPE                     *fmi3SetFMUState;
    fmi3FreeFMUStateTYPE                    *fmi3FreeFMUState;
    fmi3SerializedFMUStateSizeTYPE          *fmi3SerializedFMUStateSize;
    fmi3SerializeFMUStateTYPE               *fmi3SerializeFMUState;
    fmi3DeserializeFMUStateTYPE             *fmi3DeserializeFMUState;

    /* Getting partial derivatives */
    fmi3GetDirectionalDerivativeTYPE        *fmi3GetDirectionalDerivative;
    fmi3GetAdjointDerivativeTYPE            *fmi3GetAdjointDerivative;

    /* Entering and exiting the Configuration or Reconfiguration Mode */
    fmi3EnterConfigurationModeTYPE          *fmi3EnterConfigurationMode;
    fmi3ExitConfigurationModeTYPE           *fmi3ExitConfigurationMode;

    /* Clock related functions */
    fmi3GetIntervalDecimalTYPE              *fmi3GetIntervalDecimal;
    fmi3GetIntervalFractionTYPE             *fmi3GetIntervalFraction;
    fmi3GetShiftDecimalTYPE                 *fmi3GetShiftDecimal;
    fmi3GetShiftFractionTYPE                *fmi3GetShiftFraction;
    fmi3SetIntervalDecimalTYPE              *fmi3SetIntervalDecimal;
    fmi3SetIntervalFractionTYPE             *fmi3SetIntervalFraction;
    fmi3SetShiftDecimalTYPE                 *fmi3SetShiftDecimal;
    fmi3SetShiftFractionTYPE                *fmi3SetShiftFraction;
    fmi3EvaluateDiscreteStatesTYPE          *fmi3EvaluateDiscreteStates;
    fmi3UpdateDiscreteStatesTYPE            *fmi3UpdateDiscreteStates;

    /***************************************************
    Functions for Model Exchange
    ****************************************************/

    fmi3EnterContinuousTimeModeTYPE         *fmi3EnterContinuousTimeMode;
    fmi3CompletedIntegratorStepTYPE         *fmi3CompletedIntegratorStep;

    /* Providing independent variables and re-initialization of caching */
    fmi3SetTimeTYPE                         *fmi3SetTime;
    fmi3SetContinuousStatesTYPE             *fmi3SetContinuousStates;

    /* Evaluation of the model equations */
    fmi3GetContinuousStateDerivativesTYPE   *fmi3GetContinuousStateDerivatives;
    fmi3GetEventIndicatorsTYPE              *fmi3GetEventIndicators;
    fmi3GetContinuousStatesTYPE             *fmi3GetContinuousStates;
    fmi3GetNominalsOfContinuousStatesTYPE   *fmi3GetNominalsOfContinuousStates;
    fmi3GetNumberOfEventIndicatorsTYPE      *fmi3GetNumberOfEventIndicators;
    fmi3GetNumberOfContinuousStatesTYPE     *fmi3GetNumberOfContinuousStates;

    /***************************************************
    Functions for FMI 3.0 for Co-Simulation
    ****************************************************/

    fmi3EnterStepModeTYPE                   *fmi3EnterStepMode;
    fmi3GetOutputDerivativesTYPE            *fmi3GetOutputDerivatives;
    fmi3DoStepTYPE                          *fmi3DoStep;

    /***************************************************
    Functions for Scheduled Execution
    ****************************************************/

    fmi3ActivateModelPartitionTYPE          *fmi3ActivateModelPartition;

};


/***************************************************
Common Functions
****************************************************/

/* Inquire version numbers and setting logging status */
FMI_STATIC const char* FMI3GetVersion(FMIInstance *instance);

FMI_STATIC FMIStatus FMI3SetDebugLogging(FMIInstance *instance,
    fmi3Boolean loggingOn,
    size_t nCategories,
    const fmi3String categories[]);

/* Creation and destruction of FMU instances and setting debug status */
FMI_STATIC FMIStatus FMI3InstantiateModelExchange(
    FMIInstance *instance,
    fmi3String   instantiationToken,
    fmi3String   resourcePath,
    fmi3Boolean  visible,
    fmi3Boolean  loggingOn);

FMI_STATIC FMIStatus FMI3InstantiateCoSimulation(
    FMIInstance                   *instance,
    fmi3String                     instantiationToken,
    fmi3String                     resourcePath,
    fmi3Boolean                    visible,
    fmi3Boolean                    loggingOn,
    fmi3Boolean                    eventModeUsed,
    fmi3Boolean                    earlyReturnAllowed,
    const fmi3ValueReference       requiredIntermediateVariables[],
    size_t                         nRequiredIntermediateVariables,
    fmi3IntermediateUpdateCallback intermediateUpdate);

FMI_STATIC FMIStatus FMI3InstantiateScheduledExecution(
    FMIInstance                   *instance,
    fmi3String                     instantiationToken,
    fmi3String                     resourcePath,
    fmi3Boolean                    visible,
    fmi3Boolean                    loggingOn,
    fmi3InstanceEnvironment        instanceEnvironment,
    fmi3LogMessageCallback         logMessage,
    fmi3ClockUpdateCallback        clockUpdate,
    fmi3LockPreemptionCallback     lockPreemption,
    fmi3UnlockPreemptionCallback   unlockPreemption);

FMI_STATIC FMIStatus FMI3FreeInstance(FMIInstance *instance);

/* Enter and exit initialization mode, enter event mode, terminate and reset */
FMI_STATIC FMIStatus FMI3EnterInitializationMode(FMIInstance *instance,
    fmi3Boolean toleranceDefined,
    fmi3Float64 tolerance,
    fmi3Float64 startTime,
    fmi3Boolean stopTimeDefined,
    fmi3Float64 stopTime);

FMI_STATIC FMIStatus FMI3ExitInitializationMode(FMIInstance *instance);

FMI_STATIC FMIStatus FMI3EnterEventMode(FMIInstance *instance);

FMI_STATIC FMIStatus FMI3Terminate(FMIInstance *instance);

FMI_STATIC FMIStatus FMI3Reset(FMIInstance *instance);

/* Getting and setting variable values */
FMI_STATIC FMIStatus FMI3GetFloat32(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3Float32 values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3GetFloat64(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3Float64 values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3GetInt8(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3Int8 values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3GetUInt8(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3UInt8 values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3GetInt16(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3Int16 values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3GetUInt16(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3UInt16 values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3GetInt32(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3Int32 values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3GetUInt32(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3UInt32 values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3GetInt64(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3Int64 values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3GetUInt64(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3UInt64 values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3GetBoolean(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3Boolean values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3GetString(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3String values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3GetBinary(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    size_t sizes[],
    fmi3Binary values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3GetClock(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3Clock values[]);

FMI_STATIC FMIStatus FMI3SetFloat32(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3Float32 values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3SetFloat64(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3Float64 values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3SetInt8(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3Int8 values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3SetUInt8(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3UInt8 values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3SetInt16(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3Int16 values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3SetUInt16(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3UInt16 values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3SetInt32(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3Int32 values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3SetUInt32(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3UInt32 values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3SetInt64(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3Int64 values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3SetUInt64(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3UInt64 values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3SetBoolean(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3Boolean values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3SetString(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3String values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3SetBinary(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const size_t sizes[],
    const fmi3Binary values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3SetClock(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3Clock values[]);

/* Getting Variable Dependency Information */
FMI_STATIC FMIStatus FMI3GetNumberOfVariableDependencies(FMIInstance *instance,
    fmi3ValueReference valueReference,
    size_t* nDependencies);

FMI_STATIC FMIStatus FMI3GetVariableDependencies(FMIInstance *instance,
    fmi3ValueReference dependent,
    size_t elementIndicesOfDependent[],
    fmi3ValueReference independents[],
    size_t elementIndicesOfIndependents[],
    fmi3DependencyKind dependencyKinds[],
    size_t nDependencies);

/* Getting and setting the internal FMU state */
FMI_STATIC FMIStatus FMI3GetFMUState(FMIInstance *instance, fmi3FMUState* FMUState);

FMI_STATIC FMIStatus FMI3SetFMUState(FMIInstance *instance, fmi3FMUState  FMUState);

FMI_STATIC FMIStatus FMI3FreeFMUState(FMIInstance *instance, fmi3FMUState* FMUState);

FMI_STATIC FMIStatus FMI3SerializedFMUStateSize(FMIInstance *instance,
    fmi3FMUState FMUState,
    size_t* size);

FMI_STATIC FMIStatus FMI3SerializeFMUState(FMIInstance *instance,
    fmi3FMUState FMUState,
    fmi3Byte serializedState[],
    size_t size);

FMI_STATIC FMIStatus FMI3DeserializeFMUState(FMIInstance *instance,
    const fmi3Byte serializedState[],
    size_t size,
    fmi3FMUState* FMUState);

/* Getting partial derivatives */
FMI_STATIC FMIStatus FMI3GetDirectionalDerivative(FMIInstance *instance,
    const fmi3ValueReference unknowns[],
    size_t nUnknowns,
    const fmi3ValueReference knowns[],
    size_t nKnowns,
    const fmi3Float64 seed[],
    size_t nSeed,
    fmi3Float64 sensitivity[],
    size_t nSensitivity);

FMI_STATIC FMIStatus FMI3GetAdjointDerivative(FMIInstance *instance,
    const fmi3ValueReference unknowns[],
    size_t nUnknowns,
    const fmi3ValueReference knowns[],
    size_t nKnowns,
    const fmi3Float64 seed[],
    size_t nSeed,
    fmi3Float64 sensitivity[],
    size_t nSensitivity);

/* Entering and exiting the Configuration or Reconfiguration Mode */
FMI_STATIC FMIStatus FMI3EnterConfigurationMode(FMIInstance *instance);

FMI_STATIC FMIStatus FMI3ExitConfigurationMode(FMIInstance *instance);

FMI_STATIC FMIStatus FMI3GetIntervalDecimal(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3Float64 intervals[],
    fmi3IntervalQualifier qualifiers[]);

FMI_STATIC FMIStatus FMI3GetIntervalFraction(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3UInt64 intervalCounters[],
    fmi3UInt64 resolutions[],
    fmi3IntervalQualifier qualifiers[]);

FMI_STATIC FMIStatus FMI3GetShiftDecimal(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3Float64 shifts[]);

FMI_STATIC FMIStatus FMI3GetShiftFraction(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3UInt64 shiftCounters[],
    fmi3UInt64 resolutions[]);

FMI_STATIC FMIStatus FMI3SetIntervalDecimal(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3Float64 intervals[]);

FMI_STATIC FMIStatus FMI3SetIntervalFraction(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3UInt64 intervalCounters[],
    const fmi3UInt64 resolutions[]);

FMI_STATIC FMIStatus FMI3SetShiftDecimal(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3Float64 shifts[]);

FMI_STATIC FMIStatus FMI3SetShiftFraction(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3UInt64 shiftCounters[],
    const fmi3UInt64 resolutions[]);

FMI_STATIC FMIStatus FMI3EvaluateDiscreteStates(FMIInstance *instance);

FMI_STATIC FMIStatus FMI3UpdateDiscreteStates(FMIInstance *instance,
    fmi3Boolean* discreteStatesNeedUpdate,
    fmi3Boolean* terminateSimulation,
    fmi3Boolean* nominalsOfContinuousStatesChanged,
    fmi3Boolean* valuesOfContinuousStatesChanged,
    fmi3Boolean* nextEventTimeDefined,
    fmi3Float64* nextEventTime);

/***************************************************
Functions for Model Exchange
****************************************************/

FMI_STATIC FMIStatus FMI3EnterContinuousTimeMode(FMIInstance *instance);

FMI_STATIC FMIStatus FMI3CompletedIntegratorStep(FMIInstance *instance,
    fmi3Boolean noSetFMUStatePriorToCurrentPoint,
    fmi3Boolean* enterEventMode,
    fmi3Boolean* terminateSimulation);

/* Providing independent variables and re-initialization of caching */
FMI_STATIC FMIStatus FMI3SetTime(FMIInstance *instance, fmi3Float64 time);

FMI_STATIC FMIStatus FMI3SetContinuousStates(FMIInstance *instance,
    const fmi3Float64 continuousStates[],
    size_t nContinuousStates);

/* Evaluation of the model equations */
FMI_STATIC FMIStatus FMI3GetContinuousStateDerivatives(FMIInstance *instance,
    fmi3Float64 derivatives[],
    size_t nContinuousStates);

FMI_STATIC FMIStatus FMI3GetEventIndicators(FMIInstance *instance,
    fmi3Float64 eventIndicators[],
    size_t nEventIndicators);

FMI_STATIC FMIStatus FMI3GetContinuousStates(FMIInstance *instance,
    fmi3Float64 continuousStates[],
    size_t nContinuousStates);

FMI_STATIC FMIStatus FMI3GetNominalsOfContinuousStates(FMIInstance *instance,
    fmi3Float64 nominals[],
    size_t nContinuousStates);

FMI_STATIC FMIStatus FMI3GetNumberOfEventIndicators(FMIInstance *instance,
    size_t* nEventIndicators);

FMI_STATIC FMIStatus FMI3GetNumberOfContinuousStates(FMIInstance *instance,
    size_t* nContinuousStates);

/***************************************************
Functions for Co-Simulation
****************************************************/

/* Simulating the FMU */

FMI_STATIC FMIStatus FMI3EnterStepMode(FMIInstance *instance);

FMI_STATIC FMIStatus FMI3GetOutputDerivatives(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3Int32 orders[],
    fmi3Float64 values[],
    size_t nValues);

FMI_STATIC FMIStatus FMI3DoStep(FMIInstance *instance,
    fmi3Float64 currentCommunicationPoint,
    fmi3Float64 communicationStepSize,
    fmi3Boolean noSetFMUStatePriorToCurrentPoint,
    fmi3Boolean* eventHandlingNeeded,
    fmi3Boolean* terminateSimulation,
    fmi3Boolean* earlyReturn,
    fmi3Float64* lastSuccessfulTime);

/***************************************************
Functions for Scheduled Execution
****************************************************/

FMI_STATIC FMIStatus FMI3ActivateModelPartition(FMIInstance *instance,
    fmi3ValueReference clockReference,
    fmi3Float64 activationTime);

#ifdef __cplusplus
}  /* end of extern "C" { */
#endif


static void cb_logMessage3(fmi3InstanceEnvironment instanceEnvironment,
                           fmi3Status status,
                           fmi3String category,
                           fmi3String message) {

    if (!instanceEnvironment) return;

    FMIInstance *instance = instanceEnvironment;

    if (!instance->logMessage) return;

    instance->logMessage(instance, (FMIStatus)status, category, message);
}

#if defined(FMI3_FUNCTION_PREFIX)
#define LOAD_SYMBOL(f) \
do { \
    instance->fmi3Functions->fmi3 ## f = fmi3 ## f; \
} while (0)
#elif defined(_WIN32)
#define LOAD_SYMBOL(f) \
do { \
    instance->fmi3Functions->fmi3 ## f = (fmi3 ## f ## TYPE*)GetProcAddress(instance->libraryHandle, "fmi3" #f); \
    if (!instance->fmi3Functions->fmi3 ## f) { \
        instance->logMessage(instance, FMIFatal, "fatal", "Symbol fmi3" #f " is missing in shared library."); \
        return fmi3Fatal; \
    } \
} while (0)
#else
#define LOAD_SYMBOL(f) \
do { \
    instance->fmi3Functions->fmi3 ## f = (fmi3 ## f ## TYPE*)dlsym(instance->libraryHandle, "fmi3" #f); \
    if (!instance->fmi3Functions->fmi3 ## f) { \
        instance->logMessage(instance, FMIFatal, "fatal", "Symbol fmi3" #f " is missing in shared library."); \
        return FMIFatal; \
    } \
} while (0)
#endif

#define CALL(f) \
do { \
    const FMIStatus status = (FMIStatus)instance->fmi3Functions->fmi3 ## f (instance->component); \
    if (instance->logFunctionCall) { \
        instance->logFunctionCall(instance, status, "fmi3" #f "()"); \
    } \
    instance->status = status > instance->status ? status : instance->status; \
    return status; \
} while (0)

#define CALL_ARGS(f, m, ...) \
do { \
    const FMIStatus status = (FMIStatus)instance->fmi3Functions-> fmi3 ## f (instance->component, __VA_ARGS__); \
    if (instance->logFunctionCall) { \
        FMIClearLogMessageBuffer(instance); \
        FMIAppendToLogMessageBuffer(instance, "fmi3" #f "(" m ")", __VA_ARGS__); \
        instance->logFunctionCall(instance, status, instance->logMessageBuffer); \
    } \
    instance->status = status > instance->status ? status : instance->status; \
    return status; \
} while (0)

#define CALL_ARRAY(s, t) \
do { \
    const FMIStatus status = (FMIStatus)instance->fmi3Functions->fmi3 ## s ## t(instance->component, valueReferences, nValueReferences, values, nValues); \
    if (instance->logFunctionCall) { \
        FMIClearLogMessageBuffer(instance); \
        FMIAppendToLogMessageBuffer(instance, "fmi3" #s #t "(valueReferences={"); \
        FMIAppendArrayToLogMessageBuffer(instance, valueReferences, nValueReferences, NULL, FMIValueReferenceType); \
        FMIAppendToLogMessageBuffer(instance, "}, nValueReferences=%zu, values={", nValueReferences); \
        FMIAppendArrayToLogMessageBuffer(instance, values, nValues, NULL, FMI ## t ## Type); \
        FMIAppendToLogMessageBuffer(instance, "}, nValues=%zu)", nValues); \
        instance->logFunctionCall(instance, status, instance->logMessageBuffer); \
    } \
    instance->status = status > instance->status ? status : instance->status; \
    return status; \
} while (0)

/***************************************************
Types for Common Functions
****************************************************/

/* Inquire version numbers and setting logging status */
const char* FMI3GetVersion(FMIInstance *instance) {
    if (instance->logFunctionCall) {
        instance->logFunctionCall(instance, FMIOK, "fmi3GetVersion()");
    }
    return instance->fmi3Functions->fmi3GetVersion();
}

FMIStatus FMI3SetDebugLogging(FMIInstance *instance,
    fmi3Boolean loggingOn,
    size_t nCategories,
    const fmi3String categories[]) {

    const FMIStatus status = (FMIStatus)instance->fmi3Functions->fmi3SetDebugLogging(instance->component, loggingOn, nCategories, categories);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi3SetDebugLogging(loggingOn=%d, nCategories=%zu, categories={");
        FMIAppendArrayToLogMessageBuffer(instance, categories, nCategories, NULL, FMIStringType);
        FMIAppendToLogMessageBuffer(instance, "})");
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

static FMIStatus loadSymbols3(FMIInstance *instance) {

#if !defined(FMI_VERSION) || FMI_VERSION == 3

    instance->fmi3Functions = calloc(1, sizeof(FMI3Functions));

    if (!instance->fmi3Functions) {
        return FMIError;
    }

    instance->fmiVersion = FMIVersion3;

    /***************************************************
    Common Functions
    ****************************************************/

    /* Inquire version numbers and set debug logging */
    LOAD_SYMBOL(GetVersion);
    LOAD_SYMBOL(SetDebugLogging);

    /* Creation and destruction of FMU instances */
    LOAD_SYMBOL(InstantiateModelExchange);
    LOAD_SYMBOL(InstantiateCoSimulation);
    LOAD_SYMBOL(InstantiateScheduledExecution);
    LOAD_SYMBOL(FreeInstance);

    /* Enter and exit initialization mode, terminate and reset */
    LOAD_SYMBOL(EnterInitializationMode);
    LOAD_SYMBOL(ExitInitializationMode);
    LOAD_SYMBOL(EnterEventMode);
    LOAD_SYMBOL(Terminate);
    LOAD_SYMBOL(Reset);

    /* Getting and setting variable values */
    LOAD_SYMBOL(GetFloat32);
    LOAD_SYMBOL(GetFloat64);
    LOAD_SYMBOL(GetInt8);
    LOAD_SYMBOL(GetUInt8);
    LOAD_SYMBOL(GetInt16);
    LOAD_SYMBOL(GetUInt16);
    LOAD_SYMBOL(GetInt32);
    LOAD_SYMBOL(GetUInt32);
    LOAD_SYMBOL(GetInt64);
    LOAD_SYMBOL(GetUInt64);
    LOAD_SYMBOL(GetBoolean);
    LOAD_SYMBOL(GetString);
    LOAD_SYMBOL(GetBinary);
    LOAD_SYMBOL(GetClock);
    LOAD_SYMBOL(SetFloat32);
    LOAD_SYMBOL(SetFloat64);
    LOAD_SYMBOL(SetInt8);
    LOAD_SYMBOL(SetUInt8);
    LOAD_SYMBOL(SetInt16);
    LOAD_SYMBOL(SetUInt16);
    LOAD_SYMBOL(SetInt32);
    LOAD_SYMBOL(SetUInt32);
    LOAD_SYMBOL(SetInt64);
    LOAD_SYMBOL(SetUInt64);
    LOAD_SYMBOL(SetBoolean);
    LOAD_SYMBOL(SetString);
    LOAD_SYMBOL(SetBinary);
    LOAD_SYMBOL(SetClock);

    /* Getting Variable Dependency Information */
    LOAD_SYMBOL(GetNumberOfVariableDependencies);
    LOAD_SYMBOL(GetVariableDependencies);

    /* Getting and setting the internal FMU state */
    LOAD_SYMBOL(GetFMUState);
    LOAD_SYMBOL(SetFMUState);
    LOAD_SYMBOL(FreeFMUState);
    LOAD_SYMBOL(SerializedFMUStateSize);
    LOAD_SYMBOL(SerializeFMUState);
    LOAD_SYMBOL(DeserializeFMUState);

    /* Getting partial derivatives */
    LOAD_SYMBOL(GetDirectionalDerivative);
    LOAD_SYMBOL(GetAdjointDerivative);

    /* Entering and exiting the Configuration or Reconfiguration Mode */
    LOAD_SYMBOL(EnterConfigurationMode);
    LOAD_SYMBOL(ExitConfigurationMode);

    /* Clock related functions */
    LOAD_SYMBOL(GetIntervalDecimal);
    LOAD_SYMBOL(GetIntervalFraction);
    LOAD_SYMBOL(GetShiftDecimal);
    LOAD_SYMBOL(GetShiftFraction);
    LOAD_SYMBOL(SetIntervalDecimal);
    LOAD_SYMBOL(SetIntervalFraction);
    LOAD_SYMBOL(EvaluateDiscreteStates);
    LOAD_SYMBOL(UpdateDiscreteStates);

    /***************************************************
    Functions for Model Exchange
    ****************************************************/

    LOAD_SYMBOL(EnterContinuousTimeMode);
    LOAD_SYMBOL(CompletedIntegratorStep);

    /* Providing independent variables and re-initialization of caching */
    LOAD_SYMBOL(SetTime);
    LOAD_SYMBOL(SetContinuousStates);

    /* Evaluation of the model equations */
    LOAD_SYMBOL(GetContinuousStateDerivatives);
    LOAD_SYMBOL(GetEventIndicators);
    LOAD_SYMBOL(GetContinuousStates);
    LOAD_SYMBOL(GetNominalsOfContinuousStates);
    LOAD_SYMBOL(GetNumberOfEventIndicators);
    LOAD_SYMBOL(GetNumberOfContinuousStates);

    /***************************************************
    Functions for Co-Simulation
    ****************************************************/

    /* Simulating the FMU */
    LOAD_SYMBOL(EnterStepMode);
    LOAD_SYMBOL(GetOutputDerivatives);
    LOAD_SYMBOL(DoStep);
    LOAD_SYMBOL(ActivateModelPartition);

    instance->state = FMI2StartAndEndState;

    return FMIOK;

#else

    return FMIError;

#endif
}

/* Creation and destruction of FMU instances and setting debug status */
FMIStatus FMI3InstantiateModelExchange(
    FMIInstance *instance,
    fmi3String   instantiationToken,
    fmi3String   resourcePath,
    fmi3Boolean  visible,
    fmi3Boolean  loggingOn) {

    const FMIStatus status = loadSymbols3(instance);

    fmi3LogMessageCallback logMessage = instance->logMessage ? cb_logMessage3 : NULL;

    instance->component = instance->fmi3Functions->fmi3InstantiateModelExchange(instance->name, instantiationToken, resourcePath, visible, loggingOn, instance, logMessage);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance,
            "fmi3InstantiateModelExchange("
            "instanceName=\"%s\", "
            "instantiationToken=\"%s\", "
            "resourcePath=\"%s\", "
            "visible=%d, "
            "loggingOn=%d, "
            "instanceEnvironment=0x%p, "
            "logMessage=0x%p)",
            instance->name,
            instantiationToken,
            resourcePath,
            visible,
            loggingOn,
            instance,
            logMessage);
        instance->logFunctionCall(instance, instance->component ? FMIOK : FMIError, instance->logMessageBuffer);
    }

    if (!instance->component) {
        return FMIError;
    }

    instance->interfaceType = FMIModelExchange;
    instance->state = FMI2InstantiatedState;

    return status;
}

FMIStatus FMI3InstantiateCoSimulation(
    FMIInstance                   *instance,
    fmi3String                     instantiationToken,
    fmi3String                     resourcePath,
    fmi3Boolean                    visible,
    fmi3Boolean                    loggingOn,
    fmi3Boolean                    eventModeUsed,
    fmi3Boolean                    earlyReturnAllowed,
    const fmi3ValueReference       requiredIntermediateVariables[],
    size_t                         nRequiredIntermediateVariables,
    fmi3IntermediateUpdateCallback intermediateUpdate) {

    if (loadSymbols3(instance) != FMIOK) {
        return FMIFatal;
    }

    fmi3LogMessageCallback logMessage = instance->logMessage ? cb_logMessage3 : NULL;

    instance->component = instance->fmi3Functions->fmi3InstantiateCoSimulation(
        instance->name,
        instantiationToken,
        resourcePath,
        visible,
        loggingOn,
        eventModeUsed,
        earlyReturnAllowed,
        requiredIntermediateVariables,
        nRequiredIntermediateVariables,
        instance,
        logMessage,
        intermediateUpdate);

    instance->fmi3Functions->eventModeUsed = eventModeUsed;

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance,
            "fmi3InstantiateCoSimulation("
            "instanceName=\"%s\", "
            "instantiationToken=\"%s\", "
            "resourcePath=\"%s\", "
            "visible=%d, "
            "loggingOn=%d, "
            "eventModeUsed=%d, "
            "earlyReturnAllowed=%d, "
            "requiredIntermediateVariables=0x%p, "
            "nRequiredIntermediateVariables=%zu, "
            "instanceEnvironment=0x%p, "
            "logMessage=0x%p, "
            "intermediateUpdate=0x%p)",
            instance->name,
            instantiationToken,
            resourcePath,
            visible,
            loggingOn,
            eventModeUsed,
            earlyReturnAllowed,
            requiredIntermediateVariables,
            nRequiredIntermediateVariables,
            instance,
            logMessage,
            intermediateUpdate);
        instance->logFunctionCall(instance, instance->component ? FMIOK : FMIError, instance->logMessageBuffer);
    }

    if (!instance->component) {
        return FMIError;
    }

    instance->interfaceType = FMICoSimulation;
    instance->state = FMI2InstantiatedState;

    return FMIOK;
}

FMIStatus FMI3InstantiateScheduledExecution(
    FMIInstance                   *instance,
    fmi3String                     instantiationToken,
    fmi3String                     resourcePath,
    fmi3Boolean                    visible,
    fmi3Boolean                    loggingOn,
    fmi3InstanceEnvironment        instanceEnvironment,
    fmi3LogMessageCallback         logMessage,
    fmi3ClockUpdateCallback        clockUpdate,
    fmi3LockPreemptionCallback     lockPreemption,
    fmi3UnlockPreemptionCallback   unlockPreemption) {

    if (loadSymbols3(instance) != FMIOK) {
        return FMIError;
    }

    fmi3LogMessageCallback _logMessage = instance->logMessage ? cb_logMessage3 : NULL;

    instance->component = instance->fmi3Functions->fmi3InstantiateScheduledExecution(
        instance->name,
        instantiationToken,
        resourcePath,
        visible,
        loggingOn,
        instance,
        _logMessage,
        clockUpdate,
        lockPreemption,
        unlockPreemption);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance,
            "fmi3InstantiateScheduledExecution("
            "instanceName=\"%s\", "
            "instantiationToken=\"%s\", "
            "resourcePath=\"%s\", "
            "visible=%d, "
            "loggingOn=%d, "
            "instanceEnvironment=0x%p, "
            "logMessage=0x%p, "
            "clockUpdate=0x%p, "
            "lockPreemption=0x%p, "
            "unlockPreemption=0x%p)",
            instance->name,
            instantiationToken,
            resourcePath,
            visible,
            loggingOn,
            instance,
            _logMessage,
            clockUpdate,
            lockPreemption,
            unlockPreemption);
        instance->logFunctionCall(instance, instance->component ? FMIOK : FMIError, instance->logMessageBuffer);
    }

    if (!instance->component) {
        return FMIError;
    }

    instance->interfaceType = FMIScheduledExecution;
    instance->state = FMI2InstantiatedState;

    return FMIOK;
}

FMIStatus FMI3FreeInstance(FMIInstance *instance) {

    if (!instance) {
        return FMIError;
    }

    instance->fmi3Functions->fmi3FreeInstance(instance->component);

    instance->component = NULL;

    if (instance->logFunctionCall) {
        instance->logFunctionCall(instance, FMIOK, "fmi3FreeInstance()");
    }

    return FMIOK;
}

/* Enter and exit initialization mode, enter event mode, terminate and reset */
FMIStatus FMI3EnterInitializationMode(FMIInstance *instance,
    fmi3Boolean toleranceDefined,
    fmi3Float64 tolerance,
    fmi3Float64 startTime,
    fmi3Boolean stopTimeDefined,
    fmi3Float64 stopTime) {

    instance->state = FMI2InitializationModeState;

    CALL_ARGS(EnterInitializationMode,
        "fmi3EnterInitializationMode(toleranceDefined=%d, tolerance=%.16g, startTime=%.16g, stopTimeDefined=%d, stopTime=%.16g)",
        toleranceDefined, tolerance, startTime, stopTimeDefined, stopTime);
}

FMIStatus FMI3ExitInitializationMode(FMIInstance *instance) {

    if (instance->interfaceType == FMIModelExchange || (instance->fmiVersion == FMIVersion3 && instance->interfaceType == FMICoSimulation && instance->fmi3Functions->eventModeUsed)) {
        instance->state = FMI2EventModeState;
    } else {
        instance->state = FMI2StepCompleteState;
    }

    CALL(ExitInitializationMode);
}

FMIStatus FMI3EnterEventMode(FMIInstance *instance) {
    instance->state = FMI2EventModeState;
    CALL(EnterEventMode);
}

FMIStatus FMI3Terminate(FMIInstance *instance) {
    instance->state = FMI2TerminatedState;
    CALL(Terminate);
}

FMIStatus FMI3Reset(FMIInstance *instance) {
    instance->state = FMI2InstantiatedState;
    CALL(Reset);
}

/* Getting and setting variable values */
FMIStatus FMI3GetFloat32(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3Float32 values[],
    size_t nValues) {
    CALL_ARRAY(Get, Float32);
}

FMIStatus FMI3GetFloat64(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3Float64 values[],
    size_t nValues) {
    CALL_ARRAY(Get, Float64);
}

FMIStatus FMI3GetInt8(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3Int8 values[],
    size_t nValues) {
    CALL_ARRAY(Get, Int8);
}

FMIStatus FMI3GetUInt8(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3UInt8 values[],
    size_t nValues) {
    CALL_ARRAY(Get, UInt8);
}

FMIStatus FMI3GetInt16(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3Int16 values[],
    size_t nValues) {
    CALL_ARRAY(Get, Int16);
}

FMIStatus FMI3GetUInt16(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3UInt16 values[],
    size_t nValues) {
    CALL_ARRAY(Get, UInt16);
}

FMIStatus FMI3GetInt32(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3Int32 values[],
    size_t nValues) {
    CALL_ARRAY(Get, Int32);
}

FMIStatus FMI3GetUInt32(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3UInt32 values[],
    size_t nValues) {
    CALL_ARRAY(Get, UInt32);
}

FMIStatus FMI3GetInt64(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3Int64 values[],
    size_t nValues) {
    CALL_ARRAY(Get, Int64);
}

FMIStatus FMI3GetUInt64(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3UInt64 values[],
    size_t nValues) {
    CALL_ARRAY(Get, UInt64);
}

FMIStatus FMI3GetBoolean(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3Boolean values[],
    size_t nValues) {
    CALL_ARRAY(Get, Boolean);
}

FMIStatus FMI3GetString(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3String values[],
    size_t nValues) {
    CALL_ARRAY(Get, String);
}

FMIStatus FMI3GetBinary(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    size_t sizes[],
    fmi3Binary values[],
    size_t nValues) {

    const FMIStatus status = (FMIStatus)instance->fmi3Functions->fmi3GetBinary(instance->component, valueReferences, nValueReferences, sizes, values, nValues);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi3GetBinary(valueReferences={");
        FMIAppendArrayToLogMessageBuffer(instance, valueReferences, nValueReferences, NULL, FMIValueReferenceType);
        FMIAppendToLogMessageBuffer(instance, "}, nValueReferences=%zu, sizes={", nValueReferences);
        FMIAppendArrayToLogMessageBuffer(instance, sizes, nValueReferences, NULL, FMISizeTType);
        FMIAppendToLogMessageBuffer(instance, "}, values={");
        FMIAppendArrayToLogMessageBuffer(instance, values, nValues, sizes, FMIBinaryType);
        FMIAppendToLogMessageBuffer(instance, "}, nValues=%zu)", nValues);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI3GetClock(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3Clock values[]) {

    const FMIStatus status = (FMIStatus)instance->fmi3Functions->fmi3GetClock(instance->component, valueReferences, nValueReferences, values);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi3GetClock(valueReferences={");
        FMIAppendArrayToLogMessageBuffer(instance, valueReferences, nValueReferences, NULL, FMIValueReferenceType);
        FMIAppendToLogMessageBuffer(instance, "}, nValueReferences=%zu, values={", nValueReferences);
        FMIAppendArrayToLogMessageBuffer(instance, values, nValueReferences, NULL, FMIClockType);
        FMIAppendToLogMessageBuffer(instance, "})");
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI3SetFloat32(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3Float32 values[],
    size_t nValues)    {
    CALL_ARRAY(Set, Float32);
}

FMIStatus FMI3SetFloat64(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3Float64 values[],
    size_t nValues) {
    CALL_ARRAY(Set, Float64);
}

FMIStatus FMI3SetInt8(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3Int8 values[],
    size_t nValues) {
    CALL_ARRAY(Set, Int8);
}

FMIStatus FMI3SetUInt8(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3UInt8 values[],
    size_t nValues) {
    CALL_ARRAY(Set, UInt8);
}

FMIStatus FMI3SetInt16(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3Int16 values[],
    size_t nValues) {
    CALL_ARRAY(Set, Int16);
}

FMIStatus FMI3SetUInt16(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3UInt16 values[],
    size_t nValues) {
    CALL_ARRAY(Set, UInt16);
}

FMIStatus FMI3SetInt32(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3Int32 values[],
    size_t nValues) {
    CALL_ARRAY(Set, Int32);
}

FMIStatus FMI3SetUInt32(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3UInt32 values[],
    size_t nValues) {
    CALL_ARRAY(Set, UInt32);
}

FMIStatus FMI3SetInt64(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3Int64 values[],
    size_t nValues) {
    CALL_ARRAY(Set, Int64);
}

FMIStatus FMI3SetUInt64(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3UInt64 values[],
    size_t nValues) {
    CALL_ARRAY(Set, UInt64);
}

FMIStatus FMI3SetBoolean(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3Boolean values[],
    size_t nValues) {
    CALL_ARRAY(Set, Boolean);
}

FMIStatus FMI3SetString(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3String values[],
    size_t nValues) {
    CALL_ARRAY(Set, String);
}

FMIStatus FMI3SetBinary(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const size_t sizes[],
    const fmi3Binary values[],
    size_t nValues) {

    const FMIStatus status = (FMIStatus)instance->fmi3Functions->fmi3SetBinary(instance->component, valueReferences, nValueReferences, sizes, values, nValues);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi3SetBinary(valueReferences={");
        FMIAppendArrayToLogMessageBuffer(instance, valueReferences, nValueReferences, NULL, FMIValueReferenceType);
        FMIAppendToLogMessageBuffer(instance, "}, nValueReferences=%zu, sizes={", nValueReferences);
        FMIAppendArrayToLogMessageBuffer(instance, sizes, nValueReferences, NULL, FMISizeTType);
        FMIAppendToLogMessageBuffer(instance, "}, values={");
        FMIAppendArrayToLogMessageBuffer(instance, values, nValues, sizes, FMIBinaryType);
        FMIAppendToLogMessageBuffer(instance, "}, nValues=%zu)", nValues);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI3SetClock(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3Clock values[]) {

    const FMIStatus status = (FMIStatus)instance->fmi3Functions->fmi3SetClock(instance->component, valueReferences, nValueReferences, values);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi3SetClock(valueReferences={");
        FMIAppendArrayToLogMessageBuffer(instance, valueReferences, nValueReferences, NULL, FMIValueReferenceType);
        FMIAppendToLogMessageBuffer(instance, "}, nValueReferences=%zu, values={", nValueReferences);
        FMIAppendArrayToLogMessageBuffer(instance, values, nValueReferences, NULL, FMIClockType);
        FMIAppendToLogMessageBuffer(instance, "})");
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

/* Getting Variable Dependency Information */
FMIStatus FMI3GetNumberOfVariableDependencies(FMIInstance *instance,
    fmi3ValueReference valueReference,
    size_t* nDependencies) {
    CALL_ARGS(GetNumberOfVariableDependencies, "valueReference=%u, nDependencies=0x%p", valueReference, nDependencies);
}

FMIStatus FMI3GetVariableDependencies(FMIInstance *instance,
    fmi3ValueReference dependent,
    size_t elementIndicesOfDependent[],
    fmi3ValueReference independents[],
    size_t elementIndicesOfIndependents[],
    fmi3DependencyKind dependencyKinds[],
    size_t nDependencies) {
    CALL_ARGS(GetVariableDependencies, "dependent=%u, elementIndicesOfDependent=0x%p, independents=0x%p, elementIndicesOfIndependents=0x%p, dependencyKinds=0x%p, nDependencies=%zu",
        dependent, elementIndicesOfDependent, independents, elementIndicesOfIndependents, dependencyKinds, nDependencies);
}

/* Getting and setting the internal FMU state */
FMIStatus FMI3GetFMUState(FMIInstance *instance, fmi3FMUState* FMUState) {
    CALL_ARGS(GetFMUState, "FMUState=0x%p", FMUState);
}

FMIStatus FMI3SetFMUState(FMIInstance *instance, fmi3FMUState  FMUState) {
    CALL_ARGS(SetFMUState, "FMUState=0x%p", FMUState);
}

FMIStatus FMI3FreeFMUState(FMIInstance *instance, fmi3FMUState* FMUState) {
    CALL_ARGS(FreeFMUState, "FMUState=0x%p", FMUState);
}


FMIStatus FMI3SerializedFMUStateSize(FMIInstance *instance,
    fmi3FMUState  FMUState,
    size_t* size) {

    const FMIStatus status = (FMIStatus)instance->fmi3Functions->fmi3SerializedFMUStateSize(instance->component, FMUState, size);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi3SerializedFMUStateSize(FMUState=0x%p, size=%zu)", FMUState, *size);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI3SerializeFMUState(FMIInstance *instance,
    fmi3FMUState  FMUState,
    fmi3Byte serializedState[],
    size_t size) {
    CALL_ARGS(SerializeFMUState, "FMUstate=0x%p, serializedState=0x%p, size=%zu", FMUState, serializedState, size);
}

FMIStatus FMI3DeserializeFMUState(FMIInstance *instance,
    const fmi3Byte serializedState[],
    size_t size,
    fmi3FMUState* FMUState) {
    CALL_ARGS(DeserializeFMUState, "serializedState=0x%p, size=%zu, FMUState=0x%p", serializedState, size, FMUState);
}

/* Getting partial derivatives */
FMIStatus FMI3GetDirectionalDerivative(FMIInstance *instance,
    const fmi3ValueReference unknowns[],
    size_t nUnknowns,
    const fmi3ValueReference knowns[],
    size_t nKnowns,
    const fmi3Float64 seed[],
    size_t nSeed,
    fmi3Float64 sensitivity[],
    size_t nSensitivity) {
    CALL_ARGS(GetDirectionalDerivative,
        "unknowns=0x%p, nUnknowns=%zu, knowns=0x%p, nKnowns=%zu, seed=0x%p, nSeed=%zu, sensitivity=0x%p, nSensitivity=%zu",
        unknowns, nUnknowns, knowns, nKnowns, seed, nSeed, sensitivity, nSensitivity);
}

FMIStatus FMI3GetAdjointDerivative(FMIInstance *instance,
    const fmi3ValueReference unknowns[],
    size_t nUnknowns,
    const fmi3ValueReference knowns[],
    size_t nKnowns,
    const fmi3Float64 seed[],
    size_t nSeed,
    fmi3Float64 sensitivity[],
    size_t nSensitivity) {
    CALL_ARGS(GetAdjointDerivative,
        "unknowns=0x%p, nUnknowns=%zu, knowns=0x%p, nKnowns=%zu, seed=0x%p, nSeed=%zu, sensitivity=0x%p, nSensitivity=%zu",
        unknowns, nUnknowns, knowns, nKnowns, seed, nSeed, sensitivity, nSensitivity);
}

/* Entering and exiting the Configuration or Reconfiguration Mode */
FMIStatus FMI3EnterConfigurationMode(FMIInstance *instance) {
    CALL(EnterConfigurationMode);
}

FMIStatus FMI3ExitConfigurationMode(FMIInstance *instance) {
    CALL(ExitConfigurationMode);
}

/* Clock related functions */

FMI_STATIC FMIStatus FMI3GetIntervalDecimal(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3Float64 intervals[],
    fmi3IntervalQualifier qualifiers[]) {
    CALL_ARGS(GetIntervalDecimal,
        "valueReferences=0x%p, nValueReferences=%zu, intervals=0x%p, qualifiers=0x%p",
        valueReferences, nValueReferences, intervals, qualifiers);
}

FMIStatus FMI3GetIntervalFraction(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3UInt64 intervalCounters[],
    fmi3UInt64 resolutions[],
    fmi3IntervalQualifier qualifiers[]) {
    CALL_ARGS(GetIntervalFraction,
        "valueReferences=0x%p, nValueReferences=%zu, intervalCounters=0x%p, resolutions=0x%p, qualifiers=%d",
        valueReferences, nValueReferences, intervalCounters, resolutions, qualifiers);
}

FMIStatus FMI3GetShiftDecimal(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3Float64 shifts[]) {
    CALL_ARGS(GetShiftDecimal,
        "valueReferences=0x%p, nValueReferences=%zu, shifts=0x%p",
        valueReferences, nValueReferences, shifts);
}

FMIStatus FMI3GetShiftFraction(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3UInt64 shiftCounters[],
    fmi3UInt64 resolutions[]) {
    CALL_ARGS(GetShiftFraction,
        "valueReferences=0x%p, nValueReferences=%zu, shiftCounters=0x%p, resolutions=0x%p",
        valueReferences, nValueReferences, shiftCounters, resolutions);
}

FMIStatus FMI3SetIntervalDecimal(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3Float64 intervals[]) {
    CALL_ARGS(SetIntervalDecimal,
        "valueReferences=0x%p, nValueReferences=%zu, intervals=0x%p",
        valueReferences, nValueReferences, intervals);
}

FMIStatus FMI3SetIntervalFraction(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3UInt64 intervalCounters[],
    const fmi3UInt64 resolutions[]) {
    CALL_ARGS(SetIntervalFraction,
        "valueReferences=0x%p, nValueReferences=%zu, intervalCounters=0x%p, resolutions=0x%p",
        valueReferences, nValueReferences, intervalCounters, resolutions);
}

FMIStatus FMI3SetShiftDecimal(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3Float64 shifts[]) {
    CALL_ARGS(SetShiftDecimal,
        "valueReferences=0x%p, nValueReferences=%zu, shifts=0x%p",
        valueReferences, nValueReferences, shifts);
}

FMIStatus FMI3SetShiftFraction(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3UInt64 shiftCounters[],
    const fmi3UInt64 resolutions[]) {
    CALL_ARGS(SetShiftFraction,
        "valueReferences=0x%p, nValueReferences=%zu, shiftCounters=0x%p, resolutions=0x%p",
        valueReferences, nValueReferences, shiftCounters, resolutions);
}

FMIStatus FMI3EvaluateDiscreteStates(FMIInstance *instance) {
    CALL(EvaluateDiscreteStates);
}

FMIStatus FMI3UpdateDiscreteStates(FMIInstance *instance,
    fmi3Boolean* discreteStatesNeedUpdate,
    fmi3Boolean* terminateSimulation,
    fmi3Boolean* nominalsOfContinuousStatesChanged,
    fmi3Boolean* valuesOfContinuousStatesChanged,
    fmi3Boolean* nextEventTimeDefined,
    fmi3Float64* nextEventTime) {

    const FMIStatus status = (FMIStatus)instance->fmi3Functions->fmi3UpdateDiscreteStates(instance->component, discreteStatesNeedUpdate, terminateSimulation, nominalsOfContinuousStatesChanged, valuesOfContinuousStatesChanged, nextEventTimeDefined, nextEventTime);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance,
            "fmi3UpdateDiscreteStates(discreteStatesNeedUpdate=%d, terminateSimulation=%d, nominalsOfContinuousStatesChanged=%d, valuesOfContinuousStatesChanged=%d, nextEventTimeDefined=%d, nextEventTime=%.16g)",
            *discreteStatesNeedUpdate, *terminateSimulation, *nominalsOfContinuousStatesChanged, *valuesOfContinuousStatesChanged, *nextEventTimeDefined, *nextEventTime);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

/***************************************************
Types for Functions for Model Exchange
****************************************************/

FMIStatus FMI3EnterContinuousTimeMode(FMIInstance *instance) {
    instance->state = FMI2ContinuousTimeModeState;
    CALL(EnterContinuousTimeMode);
}

FMIStatus FMI3CompletedIntegratorStep(FMIInstance *instance,
    fmi3Boolean noSetFMUStatePriorToCurrentPoint,
    fmi3Boolean* enterEventMode,
    fmi3Boolean* terminateSimulation) {

    const FMIStatus status = (FMIStatus)instance->fmi3Functions->fmi3CompletedIntegratorStep(instance->component, noSetFMUStatePriorToCurrentPoint, enterEventMode, terminateSimulation);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance,
            "fmi3CompletedIntegratorStep(noSetFMUStatePriorToCurrentPoint=%d, enterEventMode=%d, terminateSimulation=%d)",
            noSetFMUStatePriorToCurrentPoint, *enterEventMode, *terminateSimulation);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

/* Providing independent variables and re-initialization of caching */
FMIStatus FMI3SetTime(FMIInstance *instance, fmi3Float64 time) {
    instance->time = time;
    CALL_ARGS(SetTime, "time=%.16g", time);
}

FMIStatus FMI3SetContinuousStates(FMIInstance *instance,
    const fmi3Float64 continuousStates[],
    size_t nContinuousStates) {

    const FMIStatus status = (FMIStatus)instance->fmi3Functions->fmi3SetContinuousStates(instance->component, continuousStates, nContinuousStates);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi3SetContinuousStates(continuousStates={");
        FMIAppendArrayToLogMessageBuffer(instance, continuousStates, nContinuousStates, NULL, FMIFloat64Type);
        FMIAppendToLogMessageBuffer(instance, "}, nContinuousStates=%zu)", nContinuousStates);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

/* Evaluation of the model equations */
FMIStatus FMI3GetContinuousStateDerivatives(FMIInstance *instance,
    fmi3Float64 derivatives[],
    size_t nContinuousStates) {

    const FMIStatus status = (FMIStatus)instance->fmi3Functions->fmi3GetContinuousStateDerivatives(instance->component, derivatives, nContinuousStates);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi3GetContinuousStateDerivatives(derivatives={");
        FMIAppendArrayToLogMessageBuffer(instance, derivatives, nContinuousStates, NULL, FMIFloat64Type);
        FMIAppendToLogMessageBuffer(instance, "}, nContinuousStates=%zu)", nContinuousStates);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI3GetEventIndicators(FMIInstance *instance,
    fmi3Float64 eventIndicators[],
    size_t nEventIndicators) {

    const FMIStatus status = (FMIStatus)instance->fmi3Functions->fmi3GetEventIndicators(instance->component, eventIndicators, nEventIndicators);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi3GetEventIndicators(eventIndicators={");
        FMIAppendArrayToLogMessageBuffer(instance, eventIndicators, nEventIndicators, NULL, FMIFloat64Type);
        FMIAppendToLogMessageBuffer(instance, "}, nEventIndicators=%zu)", nEventIndicators);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI3GetContinuousStates(FMIInstance *instance,
    fmi3Float64 continuousStates[],
    size_t nContinuousStates) {

    const FMIStatus status = (FMIStatus)instance->fmi3Functions->fmi3GetContinuousStates(instance->component, continuousStates, nContinuousStates);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi3GetContinuousStates(continuousStates={");
        FMIAppendArrayToLogMessageBuffer(instance, continuousStates, nContinuousStates, NULL, FMIFloat64Type);
        FMIAppendToLogMessageBuffer(instance, "}, nContinuousStates=%zu)", nContinuousStates);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}

FMIStatus FMI3GetNominalsOfContinuousStates(FMIInstance *instance,
    fmi3Float64 nominals[],
    size_t nContinuousStates) {

    const FMIStatus status = (FMIStatus)instance->fmi3Functions->fmi3GetNominalsOfContinuousStates(instance->component, nominals, nContinuousStates);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance, "fmi3GetNominalsOfContinuousStates(nominals={");
        FMIAppendArrayToLogMessageBuffer(instance, nominals, nContinuousStates, NULL, FMIFloat64Type);
        FMIAppendToLogMessageBuffer(instance, "}, nContinuousStates=%zu)", nContinuousStates);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    return status;
}


FMIStatus FMI3GetNumberOfEventIndicators(FMIInstance *instance,
    size_t* nEventIndicators) {
    CALL_ARGS(GetNumberOfEventIndicators, "nEventIndicators=%zu", nEventIndicators);
}

FMIStatus FMI3GetNumberOfContinuousStates(FMIInstance *instance,
    size_t* nContinuousStates) {
    CALL_ARGS(GetNumberOfContinuousStates, "nContinuousStates=%zu", nContinuousStates);
}

/***************************************************
Types for Functions for Co-Simulation
****************************************************/

/* Simulating the FMU */

FMIStatus FMI3EnterStepMode(FMIInstance *instance) {
    CALL(EnterStepMode);
}

FMIStatus FMI3GetOutputDerivatives(FMIInstance *instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3Int32 orders[],
    fmi3Float64 values[],
    size_t nValues) {
    CALL_ARGS(GetOutputDerivatives,
        "valueReferences=0x%p, nValueReferences=%zu, orders=0x%p, values=0x%p, nValues=%zu",
        valueReferences, nValueReferences, orders, values, nValues);
}

FMIStatus FMI3DoStep(FMIInstance *instance,
    fmi3Float64 currentCommunicationPoint,
    fmi3Float64 communicationStepSize,
    fmi3Boolean noSetFMUStatePriorToCurrentPoint,
    fmi3Boolean* eventEncountered,
    fmi3Boolean* terminate,
    fmi3Boolean* earlyReturn,
    fmi3Float64* lastSuccessfulTime) {

    const FMIStatus status = (FMIStatus)instance->fmi3Functions->fmi3DoStep(instance->component, currentCommunicationPoint, communicationStepSize, noSetFMUStatePriorToCurrentPoint, eventEncountered, terminate, earlyReturn, lastSuccessfulTime);

    if (instance->logFunctionCall) {
        FMIClearLogMessageBuffer(instance);
        FMIAppendToLogMessageBuffer(instance,
            "fmi3DoStep(currentCommunicationPoint=%.16g, communicationStepSize=%.16g, noSetFMUStatePriorToCurrentPoint=%d, eventEncountered=%d, terminate=%d, earlyReturn=%d, lastSuccessfulTime=%.16g)",
            currentCommunicationPoint, communicationStepSize, noSetFMUStatePriorToCurrentPoint, *eventEncountered, *terminate, *earlyReturn, *lastSuccessfulTime);
        instance->logFunctionCall(instance, status, instance->logMessageBuffer);
    }

    instance->time = *lastSuccessfulTime;

    return status;
}

FMIStatus FMI3ActivateModelPartition(FMIInstance *instance,
    fmi3ValueReference clockReference,
    fmi3Float64 activationTime) {
    CALL_ARGS(ActivateModelPartition,
        "clockReference=%u, activationTime=%.16g",
        clockReference, activationTime);
}

#undef LOAD_SYMBOL
#undef CALL
#undef CALL_ARGS
#undef CALL_ARRAY

#ifdef _WIN32

#include "shlwapi.h"
#pragma comment(lib, "shlwapi.lib")

#else

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

#endif

#ifndef S_FUNCTION_NAME
#define S_FUNCTION_NAME sfun_fmurun
#endif

#define S_FUNCTION_LEVEL 2

#include "simstruc.h"


typedef enum {

	fmiVersionParam,
	runAsKindParam,
	guidParam,
	modelIdentifierParam,
	unzipDirectoryParam,
    debugLoggingParam,
    logFMICallsParam,
    logLevelParam,
    logFileParam,
	relativeToleranceParam,
	sampleTimeParam,
	offsetTimeParam,
	nxParam,
	nzParam,
    resettableParam,
	inputPortWidthsParam,
	inputPortDirectFeedThroughParam,
	inputPortTypesParam,
	inputPortVariableVRsParam,
	outputPortWidthsParam,
	outputPortTypesParam,
	outputPortVariableVRsParam,
    numParams

} Parameter;

static size_t typeSizes[13] = {
    sizeof(real32_T),  //FMIFloat32Type,
    sizeof(real32_T),  //FMIDiscreteFloat32Type,
    sizeof(real_T),    //FMIFloat64Type,
    sizeof(real_T),    //FMIDiscreteFloat64Type,
    sizeof(int8_T),    //FMIInt8Type,
    sizeof(uint8_T),   //FMIUInt8Type,
    sizeof(int16_T),   //FMIInt16Type,
    sizeof(uint16_T),  //FMIUInt16Type,
    sizeof(int32_T),   //FMIInt32Type,
    sizeof(uint32_T),  //FMIUInt32Type,
    sizeof(int32_T),   //FMIInt64Type,
    sizeof(uint32_T),  //FMIUInt64Type,
    sizeof(boolean_T), //FMIBooleanType,
};

static char* getStringParam(SimStruct *S, int parameter, int index) {

	const mxArray *array = ssGetSFcnParam(S, parameter);

	const int m = (int)mxGetM(array);  // number of strings
	const int n = (int)mxGetN(array);  // max length
	const mxChar *data = (mxChar *)mxGetData(array);

	char *cstr = (char *)mxMalloc(n + 1);
	memset(cstr, '\0', n + 1);

	if (index >= m || n < 1 || !data) {
		return cstr;
	}

	// copy the row
	for (int j = 0; j < n; j++) {
		cstr[j] = (char)data[j * m + index];
	}

	// remove the trailing blanks
	for (int j = n - 1; j >= 0; j--) {
		if (cstr[j] != ' ') break;
		cstr[j] = '\0';
	}
	
	return cstr; // must be mxFree()'d
}

static bool isFMI1(SimStruct *S) {
	return mxGetScalar(ssGetSFcnParam(S, fmiVersionParam)) == 1.0;
}

static bool isFMI2(SimStruct *S) {
	return mxGetScalar(ssGetSFcnParam(S, fmiVersionParam)) == 2.0;
}

static bool isFMI3(SimStruct *S) {
	return mxGetScalar(ssGetSFcnParam(S, fmiVersionParam)) == 3.0;
}

static bool isME(SimStruct *S) { 
	return mxGetScalar(ssGetSFcnParam(S, runAsKindParam)) == FMIModelExchange;
}

static bool isCS(SimStruct *S) {
	return mxGetScalar(ssGetSFcnParam(S, runAsKindParam)) == FMICoSimulation;
}

static bool debugLogging(SimStruct *S) {
    return mxGetScalar(ssGetSFcnParam(S, debugLoggingParam));
}

static bool logFMICalls(SimStruct *S) {
    return mxGetScalar(ssGetSFcnParam(S, logFMICallsParam));
}

static FMIStatus logLevel(SimStruct *S) {
    return (int)mxGetScalar(ssGetSFcnParam(S, logLevelParam));
}

static double relativeTolerance(SimStruct *S) {
    return mxGetScalar(ssGetSFcnParam(S, relativeToleranceParam));
}

static double sampleTime(SimStruct *S) {
    return mxGetScalar(ssGetSFcnParam(S, sampleTimeParam));
}

static double offsetTime(SimStruct *S) {
    return mxGetScalar(ssGetSFcnParam(S, offsetTimeParam));
}

// number of continuous states
static int nx(SimStruct *S) {
    return (int)mxGetScalar(ssGetSFcnParam(S, nxParam));
}

// number of zero-crossings
static int nz(SimStruct *S) {
    return (int)mxGetScalar(ssGetSFcnParam(S, nzParam));
}

static bool resettable(SimStruct *S) {
    return mxGetScalar(ssGetSFcnParam(S, resettableParam));
}

static int inputPortWidth(SimStruct *S, int index) {
	const real_T *portWidths = (const real_T *)mxGetData(ssGetSFcnParam(S, inputPortWidthsParam));
	return (int)portWidths[index];
}

static bool inputPortDirectFeedThrough(SimStruct *S, int index) {
	const real_T *directFeedThrough = (const real_T *)mxGetData(ssGetSFcnParam(S, inputPortDirectFeedThroughParam));
	return directFeedThrough[index] != 0;
}

// number of input ports
static int nu(SimStruct *S) {
	return (int)mxGetNumberOfElements(ssGetSFcnParam(S, inputPortWidthsParam));
}

// number of input variables
static int nuv(SimStruct *S) {
	return (int)mxGetNumberOfElements(ssGetSFcnParam(S, inputPortVariableVRsParam));
}

static int startValueSize(SimStruct *S, Parameter parameter, int index) {
	const real_T *sizes = (const real_T *)mxGetData(ssGetSFcnParam(S, parameter));
	return (int)sizes[index];
}

static FMIValueReference valueReference(SimStruct *S, Parameter parameter, int index) {
	const mxArray *param = ssGetSFcnParam(S, parameter);
	const real_T realValue = ((const real_T *)mxGetData(param))[index];
	return (FMIValueReference)realValue;
}

static FMIVariableType variableType(SimStruct *S, Parameter parameter, int index) {
	const mxArray *param = ssGetSFcnParam(S, parameter);
	const real_T realValue = ((const real_T *)mxGetData(param))[index];
	const int intValue = (int)realValue;
	return (FMIVariableType)intValue;
}

static DTypeId simulinkVariableType(SimStruct *S, Parameter parameter, size_t index) {
	
	const mxArray *param = ssGetSFcnParam(S, parameter);
	const real_T realValue = ((real_T *)mxGetData(param))[index];
	const int intValue = (int)realValue;
	const FMIVariableType type = (FMIVariableType)intValue;
	
	switch (type) {
	case FMIFloat32Type:
	case FMIDiscreteFloat32Type:
		return SS_SINGLE;
	case FMIFloat64Type:
	case FMIDiscreteFloat64Type:
		return SS_DOUBLE;
	case FMIInt8Type:     
		return SS_INT8;
	case FMIUInt8Type:    
		return SS_UINT8;
	case FMIInt16Type:    
		return SS_INT16;
	case FMIUInt16Type:   
		return SS_UINT16;
	case FMIInt32Type:    
		return SS_INT32;
	case FMIUInt32Type:   
		return SS_UINT32;
	case FMIInt64Type:    
		return SS_INT32;
	case FMIUInt64Type:   
		return SS_UINT32;
	case FMIBooleanType:  
		return SS_BOOLEAN;
	default:              
		return -1; // error
	}
}

static real_T scalarValue(SimStruct *S, Parameter parameter, size_t index) {
	const mxArray *param = ssGetSFcnParam(S, parameter);
	return ((real_T *)mxGetData(param))[index];
}

static int outputPortWidth(SimStruct *S, size_t index) {
	const real_T *portWidths = (const real_T *)mxGetData(ssGetSFcnParam(S, outputPortWidthsParam));
	return (int)portWidths[index];
}

static int ny(SimStruct *S) { 
	return (int)mxGetNumberOfElements(ssGetSFcnParam(S, outputPortWidthsParam));
}

static bool initialized(SimStruct* S) {

	void** p = ssGetPWork(S);

	FMIInstance* instance = (FMIInstance*)p[0];

	switch (instance->state) {
	case FMI2EventModeState:
	case FMI2ContinuousTimeModeState:
	case FMI2StepCompleteState:
		return true;
	default: 
		return false;
	}
}

static void logCall(SimStruct *S, const char* message) {

    FILE *logfile = NULL;

	void **p = ssGetPWork(S);

    if (p) {
        logfile = (FILE *)p[1];
    }

    if (logfile) {
        fputs(message, logfile);
        fputs("\n", logfile);
        fflush(logfile);
    } else {
        ssPrintf(message);
        ssPrintf("\n");
    }
}

static void appendStatus(FMIStatus status, char *message, size_t size) {
	
	const char *ret;

	switch (status) {
	case FMIOK:
		ret = " -> OK";
		break;
	case FMIWarning:
		ret = " -> Warning";
		break;
	case FMIDiscard:
		ret = " -> Discard";
		break;
	case FMIError:
		ret = " -> Error";
		break;
	case FMIFatal:
		ret = " -> Fatal";
		break;
	case FMIPending:
		ret = " -> Pending";
		break;
	default:
		ret = "Illegal status code";
		break;
	}

	strncat(message, ret, size);
}

static void cb_logMessage(FMIInstance *instance, FMIStatus status, const char *category, const char * message) {
	
	SimStruct *S = (SimStruct *)instance->userData;

	if (status < logLevel(S)) {
		return;
	}

	char buf[FMI_MAX_MESSAGE_LENGTH];

	size_t len = snprintf(buf, FMI_MAX_MESSAGE_LENGTH, "[%s] ", instance->name);

	strncat(buf, message, FMI_MAX_MESSAGE_LENGTH);

	appendStatus(status, buf, FMI_MAX_MESSAGE_LENGTH);
	
	logCall(S, buf);
}

static void cb_logFunctionCall(FMIInstance *instance, FMIStatus status, const char *message) {
	
	char buf[FMI_MAX_MESSAGE_LENGTH];

	snprintf(buf, FMI_MAX_MESSAGE_LENGTH, "[%s] %s", instance->name, message);

	appendStatus(status, buf, FMI_MAX_MESSAGE_LENGTH);

	SimStruct *S = (SimStruct *)instance->userData;

	logCall(S, buf);
}

#ifdef _WIN32
const char * exceptionCodeToString(DWORD exceptionCode) {
	switch (exceptionCode) {
	case STATUS_WAIT_0:						return "WAIT_0";
	case STATUS_ABANDONED_WAIT_0:			return "ABANDONED_WAIT_0";
	case STATUS_USER_APC:					return "USER_APC";
	case STATUS_TIMEOUT:					return "TIMEOUT";
	case STATUS_PENDING:					return "PENDING";
	case DBG_EXCEPTION_HANDLED:				return "EXCEPTION_HANDLED";
	case DBG_CONTINUE:						return "CONTINUE";
	case STATUS_SEGMENT_NOTIFICATION:		return "SEGMENT_NOTIFICATION";
	case STATUS_FATAL_APP_EXIT:				return "FATAL_APP_EXIT";
	case DBG_TERMINATE_THREAD:				return "TERMINATE_THREAD";
	case DBG_TERMINATE_PROCESS:				return "TERMINATE_PROCESS";
	case DBG_CONTROL_C:						return "CONTROL_C";
	case DBG_PRINTEXCEPTION_C:				return "PRINTEXCEPTION_C";
	case DBG_RIPEXCEPTION:					return "RIPEXCEPTION";
	case DBG_CONTROL_BREAK:					return "CONTROL_BREAK";
	case DBG_COMMAND_EXCEPTION:				return "COMMAND_EXCEPTION";
	case STATUS_GUARD_PAGE_VIOLATION:		return "GUARD_PAGE_VIOLATION";
	case STATUS_DATATYPE_MISALIGNMENT:		return "DATATYPE_MISALIGNMENT";
	case STATUS_BREAKPOINT:					return "BREAKPOINT";
	case STATUS_SINGLE_STEP:				return "SINGLE_STEP";
	case STATUS_LONGJUMP:					return "LONGJUMP";
	case STATUS_UNWIND_CONSOLIDATE:			return "UNWIND_CONSOLIDATE";
	case DBG_EXCEPTION_NOT_HANDLED:			return "EXCEPTION_NOT_HANDLED";
	case STATUS_ACCESS_VIOLATION:			return "ACCESS_VIOLATION";
	case STATUS_IN_PAGE_ERROR:				return "IN_PAGE_ERROR";
	case STATUS_INVALID_HANDLE:				return "INVALID_HANDLE";
	case STATUS_INVALID_PARAMETER:			return "INVALID_PARAMETER";
	case STATUS_NO_MEMORY:					return "NO_MEMORY";
	case STATUS_ILLEGAL_INSTRUCTION:		return "ILLEGAL_INSTRUCTION";
	case STATUS_NONCONTINUABLE_EXCEPTION:	return "NONCONTINUABLE_EXCEPTION";
	case STATUS_INVALID_DISPOSITION:		return "INVALID_DISPOSITION";
	case STATUS_ARRAY_BOUNDS_EXCEEDED:		return "ARRAY_BOUNDS_EXCEEDED";
	case STATUS_FLOAT_DENORMAL_OPERAND:		return "FLOAT_DENORMAL_OPERAND";
	case STATUS_FLOAT_DIVIDE_BY_ZERO:		return "FLOAT_DIVIDE_BY_ZERO";
	case STATUS_FLOAT_INEXACT_RESULT:		return "FLOAT_INEXACT_RESULT";
	case STATUS_FLOAT_INVALID_OPERATION:	return "FLOAT_INVALID_OPERATION";
	case STATUS_FLOAT_OVERFLOW:				return "FLOAT_OVERFLOW";
	case STATUS_FLOAT_STACK_CHECK:			return "FLOAT_STACK_CHECK";
	case STATUS_FLOAT_UNDERFLOW:			return "FLOAT_UNDERFLOW";
	case STATUS_INTEGER_DIVIDE_BY_ZERO:		return "INTEGER_DIVIDE_BY_ZERO";
	case STATUS_INTEGER_OVERFLOW:			return "INTEGER_OVERFLOW";
	case STATUS_PRIVILEGED_INSTRUCTION:		return "PRIVILEGED_INSTRUCTION";
	case STATUS_STACK_OVERFLOW:				return "STACK_OVERFLOW";
	case STATUS_DLL_NOT_FOUND:				return "DLL_NOT_FOUND";
	case STATUS_ORDINAL_NOT_FOUND:			return "ORDINAL_NOT_FOUND";
	case STATUS_ENTRYPOINT_NOT_FOUND:		return "ENTRYPOINT_NOT_FOUND";
	case STATUS_CONTROL_C_EXIT:				return "CONTROL_C_EXIT";
	case STATUS_DLL_INIT_FAILED:			return "DLL_INIT_FAILED";
	case STATUS_FLOAT_MULTIPLE_FAULTS:		return "FLOAT_MULTIPLE_FAULTS";
	case STATUS_FLOAT_MULTIPLE_TRAPS:		return "FLOAT_MULTIPLE_TRAPS";
	case STATUS_REG_NAT_CONSUMPTION:		return "REG_NAT_CONSUMPTION";
	//case STATUS_HEAP_CORRUPTION:			return "HEAP_CORRUPTION";
	case STATUS_STACK_BUFFER_OVERRUN:		return "STACK_BUFFER_OVERRUN";
	case STATUS_INVALID_CRUNTIME_PARAMETER: return "INVALID_CRUNTIME_PARAMETER";
	case STATUS_ASSERTION_FAILURE:			return "ASSERTION_FAILURE";
	case STATUS_SXS_EARLY_DEACTIVATION:		return "SXS_EARLY_DEACTIVATION";
	case STATUS_SXS_INVALID_DEACTIVATION:	return "SXS_INVALID_DEACTIVATION";
	default:								return "UNKOWN_EXEPTION_CODE";
	}
}

#define CHECK_STATUS(s) \
	__try { \
		if (s > fmi2Warning) { \
			setErrorStatus(S, "The FMU reported an error."); \
			return; \
		} \
	} __except (EXCEPTION_EXECUTE_HANDLER) { \
		setErrorStatus(S, "The FMU crashed (exception code 0x%lx (%s)).", GetExceptionCode(), exceptionCodeToString(GetExceptionCode())); \
		return; \
	}
#else
#define CHECK_STATUS(s) if (s > fmi2Warning) { ssSetErrorStatus(S, "The FMU encountered an error."); return; }
#endif


#define CHECK_ERROR(f) f; if (ssGetErrorStatus(S)) return;


/* log mdl*() and fmi*() calls */
static void logDebug(SimStruct *S, const char* message, ...) {

    if (logFMICalls(S)) {
		
		char buf[FMI_MAX_MESSAGE_LENGTH];
		
		size_t len = snprintf(buf, FMI_MAX_MESSAGE_LENGTH, "[%s] ", ssGetPath(S));
		
		va_list args;
		va_start(args, message);
        vsnprintf(&buf[len], FMI_MAX_MESSAGE_LENGTH - len, message, args);
        va_end(args);

        logCall(S, buf);
    }
}

static void setErrorStatus(SimStruct *S, const char *message, ...) {
	va_list args;
	va_start(args, message);
	static char msg[1024];
	vsnprintf(msg, 1024, message, args);
	ssSetErrorStatus(S, msg);
	va_end(args);
}

static void setInput(SimStruct *S, bool direct, bool discrete, bool *inputEvent) {

    *inputEvent = false;

	void **p = ssGetPWork(S);

	FMIInstance *instance = (FMIInstance *)p[0];
    const char *preU = p[3];

	int iu  = 0;  // input port index
    int ipu = 0;  // previous input index

	for (int i = 0; i < nu(S); i++) {

		const int w = inputPortWidth(S, i);
        
        FMIVariableType type = variableType(S, inputPortTypesParam, i);
        const size_t typeSize = typeSizes[type];
        const bool discreteVariable = type != FMIFloat32Type && type != FMIFloat64Type;

		if (direct && !inputPortDirectFeedThrough(S, i)) {
			iu += w;
            ipu += w * typeSize;
			continue;
		}

		const void *y = ssGetInputPortSignal(S, i);

		if (isFMI1(S) || isFMI2(S)) {

			for (int j = 0; j < w; j++) {

				const FMIValueReference vr = valueReference(S, inputPortVariableVRsParam, iu);

                const char *value = &((const char *)y)[j * typeSize];
                char *preValue = &preU[ipu];

                ipu += typeSize;

                if (memcmp(value, preValue, typeSize)) {
                    if (!discreteVariable || (discreteVariable && discrete)) {
                        memcpy(preValue, value, typeSize);
                    }
                    *inputEvent |= discreteVariable;
                } else {
                    iu++;
                    continue;
                }

				// set the input
				if (isFMI1(S)) {

					switch (type) {
					case FMIRealType:
						CHECK_STATUS(FMI1SetReal(instance, &vr, 1, &((const real_T *)y)[j]));
						break;
					case FMIIntegerType:
						CHECK_STATUS(FMI1SetInteger(instance, &vr, 1, &((const int32_T *)y)[j]));
						break;
					case FMIBooleanType:
						CHECK_STATUS(FMI1SetBoolean(instance, &vr, 1, &((const boolean_T *)y)[j]));
						break;
					default:
						setErrorStatus(S, "Unsupported type id for FMI 1.0: %d", type);
						return;
					}

				} else {

                    if (type == FMIRealType || (type == FMIDiscreteRealType && discrete)) {
                        CHECK_STATUS(FMI2SetReal(instance, &vr, 1, (const fmi2Real *)value));
                    } else if (type == FMIIntegerType && discrete) {
                        CHECK_STATUS(FMI2SetInteger(instance, &vr, 1, (const int32_T *)value));
                    } else if (type == FMIBooleanType && discrete) {
                        const fmi2Boolean booleanValue = *value;
                        CHECK_STATUS(FMI2SetBoolean(instance, &vr, 1, &booleanValue));
                    } else {
                        setErrorStatus(S, "Unsupported type id for FMI 2.0: %d", type);
                        return;
                    }

				}

				iu++;
			}

		} else {

            const size_t nValues = inputPortWidth(S, i);
            const FMIValueReference vr = valueReference(S, inputPortVariableVRsParam, i);

            char *preValue = &preU[ipu];

            ipu += nValues * typeSize;

            if (memcmp(y, preValue, nValues * typeSize)) {
                if (!discreteVariable || (discreteVariable && discrete)) {
                    memcpy(preValue, y, nValues * typeSize);
                }
                *inputEvent |= discreteVariable;
            } else {
                continue;
            }

			switch (type) {
			case FMIFloat32Type:
				CHECK_STATUS(FMI3SetFloat32(instance, &vr, 1, (const real32_T *)y, nValues));
				break;
			case FMIFloat64Type:
				CHECK_STATUS(FMI3SetFloat64(instance, &vr, 1, (const real_T *)y, nValues));
				break;
			case FMIInt8Type:
				CHECK_STATUS(FMI3SetInt8(instance, &vr, 1, (const int8_T *)y, nValues));
				break;
			case FMIUInt8Type:
				CHECK_STATUS(FMI3SetUInt8(instance, &vr, 1, (const uint8_T *)y, nValues));
				break;
			case FMIInt16Type:
				CHECK_STATUS(FMI3SetInt16(instance, &vr, 1, (const int16_T *)y, nValues));
				break;
			case FMIUInt16Type:
				CHECK_STATUS(FMI3SetUInt16(instance, &vr, 1, (const uint16_T *)y, nValues));
				break;
			case FMIInt32Type:
				CHECK_STATUS(FMI3SetInt32(instance, &vr, 1, (const int32_T *)y, nValues));
				break;
			case FMIUInt32Type:
				CHECK_STATUS(FMI3SetUInt32(instance, &vr, 1, (const uint32_T *)y, nValues));
				break;
			case FMIInt64Type: {
				fmi3Int64* values = (fmi3Int64*)calloc(nValues, sizeof(fmi3Int64));
				for (int j = 0; j < nValues; j++) {
					values[j] = ((const int32_T*)y)[j];
				}
				CHECK_STATUS(FMI3SetInt64(instance, &vr, 1, values, nValues));
				free(values);
				break;
			}
			case FMIUInt64Type: {
				fmi3UInt64* values = (fmi3UInt64*)calloc(nValues, sizeof(fmi3UInt64));
				for (int j = 0; j < nValues; j++) {
					values[j] = ((const uint32_T*)y)[j];
				}
				CHECK_STATUS(FMI3SetUInt64(instance, &vr, 1, values, nValues));
				free(values);
				break;
			}
			case FMIBooleanType: {
				fmi3Boolean *values = (fmi3Boolean *)calloc(nValues, sizeof(fmi3Boolean));
				for (int j = 0; j < nValues; j++) {
					 values[j] = ((const boolean_T*)y)[j];
				}
				CHECK_STATUS(FMI3SetBoolean(instance, &vr, 1, values, nValues));
				free(values);
				break;
			}
			default:
				setErrorStatus(S, "Unsupported type id for FMI 3.0: %d", type);
				return;
			}

			iu++;
		}
	}
}

static void getOutput(SimStruct *S) {

	void **p = ssGetPWork(S);

	FMIInstance *instance = (FMIInstance *)p[0];

	int iy = 0;  // output port variable index

	for (int i = 0; i < ny(S); i++) {

		const FMIVariableType type = variableType(S, outputPortTypesParam, i);

		void *y = ssGetOutputPortSignal(S, i);

        if (isFMI1(S)) {

            for (int j = 0; j < outputPortWidth(S, i); j++) {

                const FMIValueReference vr = valueReference(S, outputPortVariableVRsParam, iy);

                switch (type) {
                case FMIRealType:
                case FMIDiscreteRealType:
                    CHECK_STATUS(FMI1GetReal(instance, &vr, 1, &((real_T *)y)[j]));
                    break;
                case FMIIntegerType:
                    CHECK_STATUS(FMI1GetInteger(instance, &vr, 1, &((int32_T *)y)[j]));
                    break;
                case FMIBooleanType:
                    CHECK_STATUS(FMI1GetBoolean(instance, &vr, 1, &((boolean_T *)y)[j]));
                    break;
                default:
                    setErrorStatus(S, "Unsupported type id for FMI 1.0: %d", type);
                    return;
                }

                iy++;
            }

        } else if (isFMI2(S)) {

            for (int j = 0; j < outputPortWidth(S, i); j++) {

                const FMIValueReference vr = valueReference(S, outputPortVariableVRsParam, iy);

                switch (type) {
                case FMIRealType:
                case FMIDiscreteRealType:
                    CHECK_STATUS(FMI2GetReal(instance, &vr, 1, &((real_T *)y)[j]));
                    break;
                case FMIIntegerType:
                    CHECK_STATUS(FMI2GetInteger(instance, &vr, 1, &((int32_T *)y)[j]));
                    break;
                case FMIBooleanType: {
                    fmi2Boolean booleanValue;
                    CHECK_STATUS(FMI2GetBoolean(instance, &vr, 1, &booleanValue));
                    ((boolean_T *)y)[j] = booleanValue;
                    break;
                }
                default:
                    setErrorStatus(S, "Unsupported type id for FMI 2.0: %d", type);
                    return;
                }
            
                iy++;
            }

		} else {

			const size_t nValues = outputPortWidth(S, i);
			const FMIValueReference vr = valueReference(S, outputPortVariableVRsParam, i);

			switch (type) {
            case FMIFloat32Type:
            case FMIDiscreteFloat32Type:
                CHECK_STATUS(FMI3GetFloat32(instance, &vr, 1, (real32_T *)y, nValues));
				break;
            case FMIFloat64Type:
            case FMIDiscreteFloat64Type:
                CHECK_STATUS(FMI3GetFloat64(instance, &vr, 1, (real_T *)y, nValues));
				break;
			case FMIInt8Type:
				CHECK_STATUS(FMI3GetInt8(instance, &vr, 1, (int8_T *)y, nValues));
				break;
			case FMIUInt8Type:
				CHECK_STATUS(FMI3GetUInt8(instance, &vr, 1, (uint8_T *)y, nValues));
				break;
			case FMIInt16Type:
				CHECK_STATUS(FMI3GetInt16(instance, &vr, 1, (int16_T *)y, nValues));
				break;
			case FMIUInt16Type:
				CHECK_STATUS(FMI3GetUInt16(instance, &vr, 1, (uint16_T *)y, nValues));
				break;
			case FMIInt32Type:
				CHECK_STATUS(FMI3GetInt32(instance, &vr, 1, (int32_T *)y, nValues));
				break;
			case FMIUInt32Type:
				CHECK_STATUS(FMI3GetUInt32(instance, &vr, 1, (uint32_T *)y, nValues));
				break;
			case FMIInt64Type: {
				fmi3Int64* values = (fmi3Int64*)calloc(nValues, sizeof(fmi3Int64));
				CHECK_STATUS(FMI3GetInt64(instance, &vr, 1, values, nValues));
				for (int j = 0; j < nValues; j++) {
					((int32_T*)y)[j] = values[j];
				}
				free(values);
				break;
			}
			case FMIUInt64Type: {
				fmi3UInt64* values = (fmi3UInt64*)calloc(nValues, sizeof(fmi3UInt64));
				CHECK_STATUS(FMI3GetUInt64(instance, &vr, 1, values, nValues));
				for (int j = 0; j < nValues; j++) {
					((uint32_T*)y)[j] = values[j];
				}
				free(values);
				break;
			}
			case FMIBooleanType: {
				fmi3Boolean *values = (fmi3Boolean *)calloc(nValues, sizeof(fmi3Boolean));
				CHECK_STATUS(FMI3GetBoolean(instance, &vr, 1, values, nValues));
				for (int j = 0; j < nValues; j++) {
					((boolean_T*)y)[j] = values[j];
				}
				free(values);
				break;
			}
			default:
				setErrorStatus(S, "Unsupported type id for FMI 3.0: %d", type);
				return;
			}
		}
	}
}


// all, only structural, only tunable
static void setParameters(SimStruct *S, bool structuralOnly, bool tunableOnly) {

    void **p = ssGetPWork(S);

    FMIInstance *instance = (FMIInstance *)p[0];
    FMIStatus status = FMIOK;

    int nSFcnParams = ssGetSFcnParamsCount(S);

    for (int i = numParams; i < nSFcnParams; i += 5) {

        const double            strucural = mxGetScalar(ssGetSFcnParam(S, i));
        const double            tunable   = mxGetScalar(ssGetSFcnParam(S, i + 1));
        const FMIVariableType   type      = mxGetScalar(ssGetSFcnParam(S, i + 2));
        const FMIValueReference vr        = mxGetScalar(ssGetSFcnParam(S, i + 3));

        if (structuralOnly && !strucural) continue;

        if (tunableOnly && !tunable) continue;

        if (isFMI2(S)) {

            if (instance->state == FMI2ContinuousTimeModeState) {
                CHECK_STATUS(FMI2EnterEventMode(instance));
            }

            for (int j = 0; j < 1; j++) {

                // TODO: iterate over array

                const mxArray *pa = ssGetSFcnParam(S, i + 4);

                switch (type) {
                case FMIRealType:
                case FMIDiscreteRealType: {
                    const fmi2Real value = mxGetScalar(pa);
                    CHECK_STATUS(FMI2SetReal(instance, &vr, 1, &value));
                    break;
                }
                case FMIIntegerType: {
                    const fmi2Integer value = mxGetScalar(pa);
                    CHECK_STATUS(FMI2SetInteger(instance, &vr, 1, &value));
                    break;
                }
                case FMIBooleanType: {
                    const fmi2Boolean value = mxGetScalar(pa);
                    CHECK_STATUS(FMI2SetBoolean(instance, &vr, 1, &value));
                    break;
                }
                case FMIStringType: {
                	const char    *value = getStringParam(S, i + 4, 0);
               		CHECK_STATUS(FMI2SetString(instance, &vr, 1, (const fmi2String *)&value));
                	mxFree(value);
                    break;
                }
                default:
                    setErrorStatus(S, "Unsupported type id for FMI 2.0: %d", type);
                    return;
                }
            }

            if (instance->state == FMI2EventModeState) {
                CHECK_STATUS(FMI2EnterContinuousTimeMode(instance));
            }
        }

    }

}

#define SET_VALUES(t) \
	{ \
		fmi3 ## t *values = calloc(s, sizeof(fmi3 ## t)); \
		for (int j = 0; j < nValues; j++) { \
			values[j] = (fmi3 ## t)realValues[iv + j]; \
		} \
		CHECK_STATUS(FMI3Set ## t(instance, &vr, 1, values, nValues)); \
		free(values); \
	}

static void update(SimStruct *S, bool inputEvent) {

	if (isCS(S)) {
		return;  // nothing to do
	}

	FMIInstance *instance = (FMIInstance *)ssGetPWork(S)[0];

	double time = instance->time;
	bool upcomingTimeEvent;
	double nextEventTime;

	if (isFMI1(S)) {
		upcomingTimeEvent = instance->fmi1Functions->eventInfo.upcomingTimeEvent;
		nextEventTime     = instance->fmi1Functions->eventInfo.nextEventTime;
	} else if (isFMI2(S)) {
		upcomingTimeEvent = instance->fmi2Functions->eventInfo.nextEventTimeDefined;
		nextEventTime     = instance->fmi2Functions->eventInfo.nextEventTime;
    } else {
        // TODO
		upcomingTimeEvent = false;
		nextEventTime     = 0;
    }

	// Work around for the event handling in Dymola FMUs:
	bool timeEvent = upcomingTimeEvent && time >= nextEventTime;

	if (timeEvent) {
		logDebug(S, "Time event at t=%.16g", time);
	}

	bool stepEvent = false;

	if (isFMI1(S)) {

		fmi1Boolean callEventUpdate = fmi1False;
		CHECK_STATUS(FMI1CompletedIntegratorStep(instance, &callEventUpdate))
		stepEvent = callEventUpdate;
	
	} else if (isFMI2(S)) {

		fmi2Boolean enterEventMode = fmi2False;
		fmi2Boolean terminateSimulation = fmi2False;
		CHECK_STATUS(FMI2CompletedIntegratorStep(instance, fmi2True, &enterEventMode, &terminateSimulation))
		if (terminateSimulation) {
			setErrorStatus(S, "The FMU requested to terminate the simulation.");
			return;
		}
		stepEvent = enterEventMode;
	
	} else {
	
		fmi3Boolean enterEventMode = fmi3False;
		fmi3Boolean terminateSimulation = fmi3False;
		CHECK_STATUS(FMI3CompletedIntegratorStep(instance, fmi3True, &enterEventMode, &terminateSimulation))
		if (terminateSimulation) {
			setErrorStatus(S, "The FMU requested to terminate the simulation.");
			return;
		}
		stepEvent = enterEventMode;
	}

	if (stepEvent) {
		logDebug(S, "Step event at t=%.16g\n", time);
	}

	bool stateEvent = false;

	void **p = ssGetPWork(S);

	fmi3Int32 *rootsFound = (fmi3Int32 *)p[2];

	if (nz(S) > 0) {

		real_T *prez = ssGetRWork(S);
		real_T *z = prez + nz(S);

		if (isFMI1(S)) {
			CHECK_STATUS(FMI1GetEventIndicators(instance, z, nz(S)));
		} else if (isFMI2(S)) {
			CHECK_STATUS(FMI2GetEventIndicators(instance, z, nz(S)));
		} else {
			CHECK_STATUS(FMI3GetEventIndicators(instance, z, nz(S)));
		}

		// check for state events
		for (int i = 0; i < nz(S); i++) {

			const bool rising  = (prez[i] < 0 && z[i] >= 0) || (prez[i] == 0 && z[i] > 0);
			const bool falling = (prez[i] > 0 && z[i] <= 0) || (prez[i] == 0 && z[i] < 0);

			if (rising || falling) {
				logDebug(S, "State event %s z[%d] at t=%.16g", rising ? "-\\+" : "+/-", i, instance->time);
				stateEvent = true;
				rootsFound[i] = rising ? 1 : -1;
			} else {
				rootsFound[i] = 0;
			}
		}

		// remember the current event indicators
		for (int i = 0; i < nz(S); i++) prez[i] = z[i];
	}

	if (inputEvent || timeEvent || stepEvent || stateEvent) {

		if (isFMI1(S)) {

            CHECK_ERROR(setInput(S, true, true, &inputEvent))
                
            CHECK_STATUS(FMI1EventUpdate(instance, fmi1False, &instance->fmi1Functions->eventInfo));
		
		} else if (isFMI2(S)) {

			CHECK_STATUS(FMI2EnterEventMode(instance));

            CHECK_ERROR(setInput(S, true, true, &inputEvent));

			do {
				CHECK_STATUS(FMI2NewDiscreteStates(instance, &instance->fmi2Functions->eventInfo));
				if (instance->fmi2Functions->eventInfo.terminateSimulation) {
					setErrorStatus(S, "The FMU requested to terminate the simulation.");
					return;
				}
			} while (instance->fmi2Functions->eventInfo.newDiscreteStatesNeeded);

			CHECK_STATUS(FMI2EnterContinuousTimeMode(instance));
		
		} else {

            CHECK_STATUS(FMI3EnterEventMode(instance));

            CHECK_ERROR(setInput(S, true, true, &inputEvent));

			do {
				CHECK_STATUS(FMI3UpdateDiscreteStates(instance,
					&instance->fmi3Functions->discreteStatesNeedUpdate,
					&instance->fmi3Functions->terminateSimulation,
					&instance->fmi3Functions->nominalsOfContinuousStatesChanged,
					&instance->fmi3Functions->valuesOfContinuousStatesChanged,
					&instance->fmi3Functions->nextEventTimeDefined,
					&instance->fmi3Functions->nextEventTime));

				if (instance->fmi3Functions->terminateSimulation) {
					setErrorStatus(S, "The FMU requested to terminate the simulation.");
					return;
				}
			} while (instance->fmi3Functions->discreteStatesNeedUpdate);

			CHECK_STATUS(FMI3EnterContinuousTimeMode(instance));
		}

		if (nx(S) > 0) {

			real_T *x = ssGetContStates(S);

			if (isFMI1(S)) {
				CHECK_STATUS(FMI1GetContinuousStates(instance, x, nx(S)));
			} else if (isFMI2(S)) {
				CHECK_STATUS(FMI2GetContinuousStates(instance, x, nx(S)));
			} else {
				CHECK_STATUS(FMI3GetContinuousStates(instance, x, nx(S)));
			}
		}

		if (nz(S) > 0) {

			real_T *prez = ssGetRWork(S);

			if (isFMI1(S)) {
				CHECK_STATUS(FMI1GetEventIndicators(instance, prez, nz(S)));
			} else if (isFMI2(S)) {
				CHECK_STATUS(FMI2GetEventIndicators(instance, prez, nz(S)));
			} else {
				CHECK_STATUS(FMI3GetEventIndicators(instance, prez, nz(S)));
			}
		}

		ssSetSolverNeedsReset(S);
	}
}

static bool isScalar(SimStruct *S, Parameter param) {
	const mxArray *array = ssGetSFcnParam(S, param);
	return mxIsNumeric(array) && mxGetNumberOfElements(array) == 1;
}

static bool isValidVariableType(FMIVariableType type) {
	switch (type) {
    case FMIFloat32Type:
    case FMIDiscreteFloat32Type:
    case FMIFloat64Type:
    case FMIDiscreteFloat64Type:
    case FMIInt8Type:
	case FMIUInt8Type:
	case FMIInt16Type:
	case FMIUInt16Type:
	case FMIInt32Type:
	case FMIUInt32Type:
	case FMIInt64Type:
	case FMIUInt64Type:
	case FMIBooleanType:
	case FMIStringType:
	// case FMIBinaryType:
	// case FMIClockType:
		return true;
	default:
		return false;
	}
}

static void initialize(SimStruct *S) {

    void **p = ssGetPWork(S);

    FMIInstance *instance = p[0];
    const time_T time = ssGetT(S);
    const time_T stopTime = ssGetTFinal(S);  // can be -1
    const bool toleranceDefined = relativeTolerance(S) > 0;

    if (isFMI1(S)) {

        CHECK_ERROR(setParameters(S, false, false));

        if (isCS(S)) {
            CHECK_STATUS(FMI1InitializeSlave(instance, time, stopTime > time, stopTime));
        } else {
            CHECK_STATUS(FMI1SetTime(instance, time));
            CHECK_STATUS(FMI1Initialize(instance, toleranceDefined, relativeTolerance(S), &instance->fmi1Functions->eventInfo));
            if (instance->fmi1Functions->eventInfo.terminateSimulation) {
                setErrorStatus(S, "Model requested termination at init");
                return;
            }
        }

    } else if (isFMI2(S)) {

        CHECK_ERROR(setParameters(S, false, false));
        CHECK_STATUS(FMI2SetupExperiment(instance, toleranceDefined, relativeTolerance(S), time, stopTime > time, stopTime));
        CHECK_STATUS(FMI2EnterInitializationMode(instance));

    } else {

        if (mxGetNumberOfElements(ssGetSFcnParam(S, inputPortWidthsParam)) > 0) {
            CHECK_STATUS(FMI3EnterConfigurationMode(instance));
            CHECK_ERROR(setParameters(S, true, false));
            CHECK_STATUS(FMI3ExitConfigurationMode(instance));
        }

        CHECK_ERROR(setParameters(S, false, false));

        CHECK_STATUS(FMI3EnterInitializationMode(instance, toleranceDefined, relativeTolerance(S), time, stopTime > time, stopTime));

    }

    if (isME(S)) {

		if (isFMI2(S)) {
			CHECK_STATUS(FMI2ExitInitializationMode(instance));
		} else {
			CHECK_STATUS(FMI3ExitInitializationMode(instance));
		}

        // initialize the continuous states
        real_T *x = ssGetContStates(S);

        if (nx(S) > 0) {

            if (isFMI1(S)) {
                CHECK_STATUS(FMI1GetContinuousStates(instance, x, nx(S)));
            } else if (isFMI2(S)) {
                CHECK_STATUS(FMI2GetContinuousStates(instance, x, nx(S)));
            } else {
                CHECK_STATUS(FMI3GetContinuousStates(instance, x, nx(S)));
            }
        }

        // initialize the event indicators
        if (nz(S) > 0) {

            real_T *prez = ssGetRWork(S);

            if (isFMI1(S)) {
                CHECK_STATUS(FMI1GetEventIndicators(instance, prez, nz(S)));
            } else if (isFMI2(S)) {
                CHECK_STATUS(FMI2GetEventIndicators(instance, prez, nz(S)));
            } else {
                CHECK_STATUS(FMI3GetEventIndicators(instance, prez, nz(S)));
            }

            real_T *z = prez + nz(S);

            memcpy(z, prez, nz(S) * sizeof(real_T));
        }
    }
}

#define MDL_ENABLE
static void mdlEnable(SimStruct *S) {
    
    logDebug(S, "mdlEnable()");

    void **p = ssGetPWork(S);

    FMIInstance *instance = p[0];

    if (instance) {

        if (isCS(S)) {

            if (isFMI1(S)) {
                CHECK_STATUS(FMI1ResetSlave(instance));
            } else if (isFMI2(S)) {
                CHECK_STATUS(FMI2Reset(instance));
            } else {
                CHECK_STATUS(FMI3Reset(instance));
            }

            CHECK_ERROR(initialize(S));
        }

        return;
    }

    if (nz(S) > 0) {
        p[2] = calloc(nz(S), sizeof(fmi3Int32)); // rootsFound
    }

    const char_T* instanceName = ssGetPath(S);

    const bool loggingOn = debugLogging(S);

    const char *modelIdentifier = getStringParam(S, modelIdentifierParam, 0);

#ifdef GRTFMI
    char *unzipdir = (char *)mxMalloc(MAX_PATH);
    strncpy(unzipdir, FMU_RESOURCES_DIR, MAX_PATH);
    strncat(unzipdir, "/", MAX_PATH);
    strncat(unzipdir, modelIdentifier, MAX_PATH);
#else
    char *unzipdir = getStringParam(S, unzipDirectoryParam, 0);
#endif

#ifdef _WIN32
    char libraryPath[MAX_PATH];
    strncpy(libraryPath, unzipdir, MAX_PATH);
    PathAppend(libraryPath, "binaries");
    if (isFMI1(S) || isFMI2(S)) {
#ifdef _WIN64
        PathAppend(libraryPath, "win64");
#else
        PathAppend(libraryPath, "win32");
#endif
    } else {
#ifdef _WIN64
        PathAppend(libraryPath, "x86_64-windows");
#else
        PathAppend(libraryPath, "i686-windows");
#endif
    }
    PathAppend(libraryPath, modelIdentifier);
    strncat(libraryPath, ".dll", MAX_PATH);
#elif defined(__APPLE__)
    char libraryPath[PATH_MAX];
    strncpy(libraryPath, unzipdir, PATH_MAX);
    strncat(libraryPath, "/binaries/darwin64/", PATH_MAX);
    strncat(libraryPath, modelIdentifier, PATH_MAX);
    strncat(libraryPath, ".dylib", PATH_MAX);
#else
    char libraryPath[PATH_MAX];
    strncpy(libraryPath, unzipdir, PATH_MAX);
    strncat(libraryPath, "/binaries/linux64/", PATH_MAX);
    strncat(libraryPath, modelIdentifier, PATH_MAX);
    strncat(libraryPath, ".so", PATH_MAX);
#endif

    instance = FMICreateInstance(instanceName, libraryPath, cb_logMessage, logFMICalls(S) ? cb_logFunctionCall : NULL);

    if (!instance) {
        setErrorStatus("Failed to load %s.", libraryPath);
        return;
    }

    instance->userData = S;

    p[0] = instance;

    const char *guid = getStringParam(S, guidParam, 0);

    char fmuResourceLocation[INTERNET_MAX_URL_LENGTH];

    if (isFMI3(S)) {
        strncpy(fmuResourceLocation, unzipdir, INTERNET_MAX_URL_LENGTH);
    } else {
#ifdef _WIN32
        DWORD fmuLocationLength = INTERNET_MAX_URL_LENGTH;
        if (UrlCreateFromPath(unzipdir, fmuResourceLocation, &fmuLocationLength, 0) != S_OK) {
            setErrorStatus(S, "Failed to create fmuResourceLocation.");
            return;
        }
#else
        strcpy(fmuResourceLocation, "file://");
        strcat(fmuResourceLocation, unzipdir);
#endif
    }

    if (!isFMI1(S)) {
        strcat(fmuResourceLocation, "/resources");
    }
    
    // instantiate the FMU
    if (isFMI1(S)) {

        if (isCS(S)) {
            CHECK_STATUS(FMI1InstantiateSlave(instance, modelIdentifier, guid, fmuResourceLocation, "application/x-fmu-sharedlibrary", 0, fmi1False, fmi1False, loggingOn));
        } else {
            CHECK_STATUS(FMI1InstantiateModel(instance, modelIdentifier, guid, loggingOn));
        }

    } else if (isFMI2(S)) {

        CHECK_STATUS(FMI2Instantiate(instance, fmuResourceLocation, isCS(S) ? fmi2CoSimulation : fmi2ModelExchange, guid, fmi2False, loggingOn));

    } else {

        if (isME(S)) {
            CHECK_STATUS(FMI3InstantiateModelExchange(instance, guid, fmuResourceLocation, fmi3False, loggingOn));
        } else {
            CHECK_STATUS(FMI3InstantiateCoSimulation(instance, guid, fmuResourceLocation, fmi3False, loggingOn, fmi3False, fmi3False, NULL, 0, NULL));
        }

    }

    // initialize the FMU instance
    CHECK_ERROR(initialize(S));

    // free string parameters
    mxFree((void *)modelIdentifier);
    mxFree((void *)unzipdir);
    mxFree((void *)guid);
}

#define MDL_DISABLE
static void mdlDisable(SimStruct *S) {
    logDebug(S, "mdlDisable()");
}

#define MDL_CHECK_PARAMETERS
#if defined(MDL_CHECK_PARAMETERS) && defined(MATLAB_MEX_FILE)
static void mdlCheckParameters(SimStruct *S) {

	logDebug(S, "mdlCheckParameters()");

	if (!isScalar(S, fmiVersionParam) || !(isFMI1(S) || isFMI2(S) || isFMI3(S))) {
		setErrorStatus(S, "Parameter %d (FMI version) must be one of 1.0, 2.0 or 3.0.", fmiVersionParam + 1);
		return;
	}

	if (!isScalar(S, runAsKindParam) || (!isME(S) && !isCS(S))) {
        setErrorStatus(S, "Parameter %d (run as kind) must be one of 0 (= Model Exchange) or 1 (= Co-Simulation).", runAsKindParam + 1);
        return;
    }

	if (!mxIsChar(ssGetSFcnParam(S, guidParam))) {
        setErrorStatus(S, "Parameter %d (GUID) must be a string.", guidParam + 1);
        return;
    }

	if (!mxIsChar(ssGetSFcnParam(S, modelIdentifierParam))) {
        setErrorStatus(S, "Parameter %d (model identifier) must be a string.", modelIdentifierParam + 1);
        return;
    }

	if (!mxIsChar(ssGetSFcnParam(S, unzipDirectoryParam))) {
		setErrorStatus(S, "Parameter %d (unzip directory) must be a string.", unzipDirectoryParam + 1);
		return;
	}

    if (!isScalar(S, debugLoggingParam)) {
        setErrorStatus(S, "Parameter %d (debug logging) must be a scalar.", debugLoggingParam + 1);
        return;
    }
	
	if (!isScalar(S, logFMICallsParam)) {
		setErrorStatus(S, "Parameter %d (log FMI calls) must be a scalar.", logFMICallsParam + 1);
		return;
	}

	if (!isScalar(S, logLevelParam) ||
		(logLevel(S) != 0 && logLevel(S) != 1 && logLevel(S) != 2 && logLevel(S) != 3 && logLevel(S) != 4 && logLevel(S) != 5)) {
		setErrorStatus(S, "Parameter %d (log level) must be one of 0 (= info), 1 (= warning), 2 (= discard), 3 (= error), 4 (= fatal) or 5 (= none).", logLevelParam + 1);
		return;
	}

	if (!mxIsChar(ssGetSFcnParam(S, logFileParam))) {
		setErrorStatus(S, "Parameter %d (log file) must be a string.", logFileParam + 1);
		return;
	}

	if (!isScalar(S, relativeToleranceParam)) {
		setErrorStatus(S, "Parameter %d (relative tolerance) must be a scalar.", relativeToleranceParam + 1);
		return;
	}

	if (!isScalar(S, sampleTimeParam)) {
		setErrorStatus(S, "Parameter %d (sample time) must be a scalar.", sampleTimeParam + 1);
		return;
	}

	if (!isScalar(S, offsetTimeParam)) {
		setErrorStatus(S, "Parameter %d (offset time) must be a scalar.", offsetTimeParam + 1);
		return;
	}

	if (!isScalar(S, nxParam)) {
		setErrorStatus(S, "Parameter %d (number of continuous states) must be a scalar.", nxParam + 1);
		return;
	}

	if (!isScalar(S, nzParam)) {
		setErrorStatus(S, "Parameter %d (number of event indicators) must be a scalar.", nzParam + 1);
		return;
	}

	if (!mxIsDouble(ssGetSFcnParam(S, inputPortWidthsParam))) {
		setErrorStatus(S, "Parameter %d (input port widths) must be a double array", inputPortWidthsParam + 1);
		return;
	}

	if (!mxIsDouble(ssGetSFcnParam(S, inputPortDirectFeedThroughParam))) {
		setErrorStatus(S, "Parameter %d (input port direct feed through) must be a double array", inputPortDirectFeedThroughParam + 1);
		return;
	}

	if (mxGetNumberOfElements(ssGetSFcnParam(S, inputPortDirectFeedThroughParam)) != mxGetNumberOfElements(ssGetSFcnParam(S, inputPortWidthsParam))) {
		setErrorStatus(S, "The number of elements in parameter %d (input port direct feed through) must be equal to the number of elements in parameter %d (inport port widths)", inputPortDirectFeedThroughParam + 1, inputPortWidthsParam + 1);
		return;
	}

	int nu = 0; // number of input variables

	if (isFMI1(S) || isFMI2(S)) {
		for (int i = 0; i < mxGetNumberOfElements(ssGetSFcnParam(S, inputPortWidthsParam)); i++) {
			if (inputPortWidth(S, i) < 1) {
				setErrorStatus(S, "Elements in parameter %d (input port widths) must be >= 1", inputPortWidthsParam + 1);
				return;
			}
			nu += inputPortWidth(S, i);
		}
	} else {
		nu = (int)mxGetNumberOfElements(ssGetSFcnParam(S, inputPortWidthsParam));
	}

	if (!mxIsDouble(ssGetSFcnParam(S, inputPortDirectFeedThroughParam)) || mxGetNumberOfElements(ssGetSFcnParam(S, inputPortDirectFeedThroughParam)) != mxGetNumberOfElements(ssGetSFcnParam(S, inputPortWidthsParam))) {
		setErrorStatus(S, "Parameter %d (input port direct feed through) must be a double array with the same number of elements as parameter %d (input port widths)", inputPortDirectFeedThroughParam + 1, inputPortWidthsParam + 1);
		return;
	}

	if (mxGetNumberOfElements(ssGetSFcnParam(S, inputPortDirectFeedThroughParam)) != mxGetNumberOfElements(ssGetSFcnParam(S, inputPortWidthsParam))) {
		setErrorStatus(S, "The number of elements in parameter %d (inport port direct feed through) must be equal to the number of elements in parameter %d (inport port widths)", inputPortDirectFeedThroughParam + 1, inputPortWidthsParam + 1);
		return;
	}

	if (!mxIsDouble(ssGetSFcnParam(S, inputPortTypesParam))) {
		setErrorStatus(S, "Parameter %d (input port variable types) must be a double array", inputPortTypesParam + 1);
		return;
	}

	for (int i = 0; i < mxGetNumberOfElements(ssGetSFcnParam(S, inputPortTypesParam)); i++) {
		FMIVariableType t = variableType(S, inputPortTypesParam, i);
		if (!isValidVariableType(t)) {
			setErrorStatus(S, "Elements in parameter %d (input port types) must be valid variable types", inputPortTypesParam + 1);
			return;
		}
	}

	if (mxGetNumberOfElements(ssGetSFcnParam(S, inputPortTypesParam)) != mxGetNumberOfElements(ssGetSFcnParam(S, inputPortWidthsParam))) {
		setErrorStatus(S, "The number of elements in parameter %d (inport port types) must be equal to the number of the elements in parameter %d (inport port widths)", inputPortTypesParam + 1, inputPortWidthsParam + 1);
		return;
	}

	if (!mxIsDouble(ssGetSFcnParam(S, inputPortVariableVRsParam))) {
		setErrorStatus(S, "Parameter %d (input port value references) must be a double array", inputPortVariableVRsParam + 1);
		return;
	}

	if (mxGetNumberOfElements(ssGetSFcnParam(S, inputPortVariableVRsParam)) != nu) {
		setErrorStatus(S, "The number of elements in parameter %d (input port value references) must be equal to the sum of the elements in parameter %d (inport port widths)", inputPortVariableVRsParam + 1, inputPortWidthsParam + 1);
		return;
	}

	// TODO: check VRS values!

	if (!mxIsDouble(ssGetSFcnParam(S, outputPortWidthsParam))) {
		setErrorStatus(S, "Parameter %d (output port widths) must be a double array", outputPortWidthsParam + 1);
		return;
	}

	int ny = 0; // number of output variables

	if (isFMI1(S) || isFMI2(S)) {
		for (int i = 0; i < mxGetNumberOfElements(ssGetSFcnParam(S, outputPortWidthsParam)); i++) {
			if (outputPortWidth(S, i) < 1) {
				setErrorStatus(S, "Elements in parameter %d (output port widths) must be >= 1", outputPortWidthsParam + 1);
				return;
			}
			ny += outputPortWidth(S, i);
		}
	} else {
		ny = (int)mxGetNumberOfElements(ssGetSFcnParam(S, outputPortWidthsParam));
	}

	if (!mxIsDouble(ssGetSFcnParam(S, outputPortTypesParam))) {
		setErrorStatus(S, "Parameter %d (output port types) must be a double array", outputPortTypesParam + 1);
		return;
	}

	if (mxGetNumberOfElements(ssGetSFcnParam(S, outputPortTypesParam)) != mxGetNumberOfElements(ssGetSFcnParam(S, outputPortWidthsParam))) {
		setErrorStatus(S, "The number of elements in parameter %d (output port types) must be equal to the number of the elements in parameter %d (output port widths)", outputPortTypesParam + 1, outputPortWidthsParam + 1);
		return;
	}

	for (int i = 0; i < mxGetNumberOfElements(ssGetSFcnParam(S, outputPortWidthsParam)); i++) {
		FMIVariableType t = variableType(S, outputPortTypesParam, i);
		if (!isValidVariableType(t)) {
			setErrorStatus(S, "Elements in parameter %d (output port types) must be valid variable types", outputPortTypesParam + 1);
			return;
		}
	}

	if (!mxIsDouble(ssGetSFcnParam(S, outputPortVariableVRsParam))) {
		setErrorStatus(S, "Parameter %d (output variable value references) must be a double array", outputPortVariableVRsParam + 1);
		return;
	}

	if (mxGetNumberOfElements(ssGetSFcnParam(S, outputPortVariableVRsParam)) != ny) {
		setErrorStatus(S, "The number of elements in parameter %d (output variable value references) must be equal to the sum of the elements in parameter %d (output port widths)", outputPortVariableVRsParam + 1, outputPortWidthsParam + 1);
		return;
	}

	// TODO: check VRS values!

}
#endif /* MDL_CHECK_PARAMETERS */


static int inputSize(SimStruct *S) {

    size_t s = 0;

    for (int i = 0; i < nu(S); i++) {
        const size_t nValues = inputPortWidth(S, i);
        FMIVariableType type = variableType(S, inputPortTypesParam, i);        
        s += nValues * typeSizes[type];
    }

    return s;
}


#define MDL_PROCESS_PARAMETERS   /* Change to #undef to remove function */
#if defined(MDL_PROCESS_PARAMETERS) && defined(MATLAB_MEX_FILE)
static void mdlProcessParameters(SimStruct *S) {

    logDebug(S, "mdlProcessParameters()");

    CHECK_ERROR(setParameters(S, false, true));
}
#endif /* MDL_PROCESS_PARAMETERS */


static void mdlInitializeSizes(SimStruct *S) {

	logDebug(S, "mdlInitializeSizes()");

    const int nSFcnParams = ssGetSFcnParamsCount(S);

    if ((nSFcnParams - numParams) % 5) {
        setErrorStatus(S, "Wrong number of arguments.");
        return;
    }

	ssSetNumSFcnParams(S, nSFcnParams);

    for (int i = 0; i < numParams; i++) {
        ssSetSFcnParamTunable(S, i, false);
    }

    for (int i = numParams; i < nSFcnParams; i += 5) {
        const double paramTunable = mxGetScalar(ssGetSFcnParam(S, i + 1));
        ssSetSFcnParamTunable(S, i, paramTunable != 0);
    }

#if defined(MATLAB_MEX_FILE)
	if (ssGetNumSFcnParams(S) == ssGetSFcnParamsCount(S)) {
		mdlCheckParameters(S);
		if (ssGetErrorStatus(S) != NULL) {
			return;
		}
	} else {
		return; // parameter mismatch will be reported by Simulink
	}
#endif

	ssSetNumContStates(S, isME(S) ? nx(S) : 0);
	ssSetNumDiscStates(S, 0);
    
    const int_T numInputPorts = nu(S) + (resettable(S) ? 1 : 0);

	if (!ssSetNumInputPorts(S, numInputPorts)) return;

	for (int i = 0; i < nu(S); i++) {
		ssSetInputPortWidth(S, i, inputPortWidth(S, i));
		ssSetInputPortRequiredContiguous(S, i, 1); // direct input signal access
		const DTypeId type = simulinkVariableType(S, inputPortTypesParam, i);
		if (type < 0) {
			setErrorStatus(S, "Unexpected type id for input port %d.", i);
			return;
		}
		ssSetInputPortDataType(S, i, type);
		bool dirFeed = inputPortDirectFeedThrough(S, i);
		ssSetInputPortDirectFeedThrough(S, i, dirFeed); // direct feed through
		logDebug(S, "ssSetInputPortDirectFeedThrough(port=%d, dirFeed=%d)", i, dirFeed);
	}

    if (resettable(S)) {
        ssSetInputPortWidth(S, nu(S), 1);
        ssSetInputPortRequiredContiguous(S, nu(S), true); // direct input signal access
        ssSetInputPortDataType(S, nu(S), SS_DOUBLE);
        ssSetInputPortDirectFeedThrough(S, nu(S), true);
    }

	if (!ssSetNumOutputPorts(S, ny(S))) return;

	for (int i = 0; i < ny(S); i++) {
		ssSetOutputPortWidth(S, i, outputPortWidth(S, i));
		DTypeId type = simulinkVariableType(S, outputPortTypesParam, i);
		if (type < 0) {
			setErrorStatus(S, "Unexpected type id for output port %d.", i);
			return;
		}
		ssSetOutputPortDataType(S, i, type);
	}

	ssSetNumSampleTimes(S, 1);
	ssSetNumRWork(S, 2 * nz(S) + nuv(S) + (resettable(S) ? 1 : 0)); // [pre(z), z, pre(u), pre(reset)]
    ssSetNumPWork(S, 4); // [FMU, logfile, rootsFound, preInput]
    ssSetNumModes(S, 3); // [stateEvent, timeEvent, stepEvent]
	ssSetNumNonsampledZCs(S, (isME(S)) ? nz(S) + 1 : 0);

	// specify the sim state compliance to be same as a built-in block
	//ssSetSimStateCompliance(S, USE_DEFAULT_SIM_STATE);

	ssSetOptions(S, 0);
}


static void mdlInitializeSampleTimes(SimStruct *S) {

	logDebug(S, "mdlInitializeSampleTimes()");

	ssSetSampleTime(S, 0, sampleTime(S));
	ssSetOffsetTime(S, 0, offsetTime(S));
}


#define MDL_START
#if defined(MDL_START)
static void mdlStart(SimStruct *S) {

    void **p = ssGetPWork(S);

    if (p[1]) {
        fclose((FILE *)p[1]);
        p[1] = NULL;
    }

	const char *logFile = getStringParam(S, logFileParam, 0);

	if (strlen(logFile) > 0) {
	    p[1] = fopen(logFile, "w");
	}

    mxFree((void *)logFile);
	if (nz(S) > 0) {
		p[2] = calloc(nz(S), sizeof(fmi3Int32));
	} else {
		p[2] = NULL;
	}

    const size_t s = inputSize(S);
    
    if (s > 0) {
        p[3] = malloc(s);
        memset(p[3], 0xFF, s);
    }

	logDebug(S, "mdlStart()");
}
#endif /* MDL_START */


static void mdlOutputs(SimStruct *S, int_T tid) {

    const time_T time = ssGetT(S);

	logDebug(S, "mdlOutputs(tid=%d, time=%.16g, majorTimeStep=%d)", tid, time, ssIsMajorTimeStep(S));

	void **p = ssGetPWork(S);

	FMIInstance *instance = (FMIInstance *)p[0];

    if (resettable(S)) {
        real_T *preReset = &(ssGetRWork(S)[2 * nz(S) + nuv(S)]);
        const real_T reset = *ssGetInputPortRealSignal(S, nu(S));

        if (*preReset == 0 && reset > 0) {
            
            if (isFMI1(S) && isCS(S)) {
                CHECK_STATUS(FMI1ResetSlave(instance));
            } else if (isFMI2(S)) {
                CHECK_STATUS(FMI2Reset(instance));
            } else if (isFMI3(S)) {
                CHECK_STATUS(FMI3Reset(instance));
            }

            CHECK_ERROR(initialize(S));
		}

        *preReset = reset;
    }

	if (isME(S)) {

		const real_T *x = ssGetContStates(S);

		if (isFMI1(S)) {

			CHECK_STATUS(FMI1SetTime(instance, time));

			if (nx(S) > 0) {
				CHECK_STATUS(FMI1SetContinuousStates(instance, x, nx(S)));
			}

		} else if (isFMI2(S)) {

			if (instance->state == FMI2EventModeState) {

                bool inputEvent;

				CHECK_ERROR(setInput(S, true, true, &inputEvent));

				do {
					CHECK_STATUS(FMI2NewDiscreteStates(instance, &instance->fmi2Functions->eventInfo));
					if (instance->fmi2Functions->eventInfo.terminateSimulation) {
						setErrorStatus(S, "The FMU requested to terminate the simulation.");
						return;
					}
				} while (instance->fmi2Functions->eventInfo.newDiscreteStatesNeeded);

				CHECK_STATUS(FMI2EnterContinuousTimeMode(instance));
			}

			if (instance->state != FMI2ContinuousTimeModeState) {
				CHECK_STATUS(FMI2EnterContinuousTimeMode(instance));
			}

			CHECK_STATUS(FMI2SetTime(instance, time));

			if (nx(S) > 0) {
				CHECK_STATUS(FMI2SetContinuousStates(instance, x, nx(S)));
			}

		} else {
			
			if (instance->state == FMI2EventModeState) {

                bool inputEvent;

                CHECK_ERROR(setInput(S, true, true, &inputEvent));

				do {
					CHECK_STATUS(FMI3UpdateDiscreteStates(instance,
						&instance->fmi3Functions->discreteStatesNeedUpdate,
						&instance->fmi3Functions->terminateSimulation,
						&instance->fmi3Functions->nominalsOfContinuousStatesChanged,
						&instance->fmi3Functions->valuesOfContinuousStatesChanged,
						&instance->fmi3Functions->nextEventTimeDefined,
						&instance->fmi3Functions->nextEventTime))

					if (instance->fmi3Functions->terminateSimulation) {
						setErrorStatus(S, "The FMU requested to terminate the simulation.");
						return;
					}
				} while (instance->fmi3Functions->discreteStatesNeedUpdate);

				CHECK_STATUS(FMI3EnterContinuousTimeMode(instance));
			}

			if (instance->state != FMI2ContinuousTimeModeState) {
				CHECK_STATUS(FMI3EnterContinuousTimeMode(instance));
			}

			CHECK_STATUS(FMI3SetTime(instance, time));

			if (nx(S) > 0) {
				CHECK_STATUS(FMI3SetContinuousStates(instance, x, nx(S)));
			}
		}

        bool inputEvent;
			   		
		CHECK_ERROR(setInput(S, true, false, &inputEvent));

		if (ssIsMajorTimeStep(S)) {
			CHECK_ERROR(update(S, inputEvent));
		}

	} else {

		const time_T h = time - instance->time;

		if (h > 0) {

			if (!initialized(S)) {
				if (isFMI2(S)) {
					CHECK_STATUS(FMI2ExitInitializationMode(instance));
				} else if (isFMI3(S)) {
					CHECK_STATUS(FMI3ExitInitializationMode(instance));
				}
			}

			if (isFMI1(S)) {
				CHECK_STATUS(FMI1DoStep(instance, instance->time, h, fmi1True));
			} else if (isFMI2(S)) {
				CHECK_STATUS(FMI2DoStep(instance, instance->time, h, fmi2True));
			} else {
				fmi3Boolean eventEncountered;
				fmi3Boolean terminateSimulation;
				fmi3Boolean earlyReturn;
				fmi3Float64 lastSuccessfulTime;
				CHECK_STATUS(FMI3DoStep(instance, instance->time, h, fmi3True, &eventEncountered, &terminateSimulation, &earlyReturn, &lastSuccessfulTime));
                // TODO: handle terminateSimulation == true
			}
		}
	}

	CHECK_ERROR(getOutput(S));
}


#define MDL_UPDATE
#if defined(MDL_UPDATE)
static void mdlUpdate(SimStruct *S, int_T tid) {

	void** p = ssGetPWork(S);

	FMIInstance* instance = (FMIInstance*)p[0];

	logDebug(S, "mdlUpdate(tid=%d, time=%.16g, majorTimeStep=%d)", tid, ssGetT(S), ssIsMajorTimeStep(S));

    bool inputEvent;

    CHECK_ERROR(setInput(S, false, isCS(S), &inputEvent));

    if (isME(S) && inputEvent) {
        ssSetErrorStatus(S, "Unexpected input event in mdlUpdate().");
    }
}
#endif // MDL_UPDATE


#define MDL_ZERO_CROSSINGS
#if defined(MDL_ZERO_CROSSINGS) && (defined(MATLAB_MEX_FILE) || defined(NRT))
static void mdlZeroCrossings(SimStruct *S) {

	logDebug(S, "mdlZeroCrossings(time=%.16g, majorTimeStep=%d)", ssGetT(S), ssIsMajorTimeStep(S));

	if (isME(S)) {

		bool inputEvent;

        CHECK_ERROR(setInput(S, true, false, &inputEvent));

        if (inputEvent) {
            ssSetErrorStatus(S, "Unexpected input event in mdlZeroCrossings().");
            return;
        }

		real_T *z = ssGetNonsampledZCs(S);

		void **p = ssGetPWork(S);

		FMIInstance *instance = (FMIInstance *)p[0];
		real_T nextEventTime;

		if (nz(S) > 0) {
			if (isFMI1(S)) {
				CHECK_STATUS(FMI1GetEventIndicators(instance, z, nz(S)));
			} else if (isFMI2(S)) {
				CHECK_STATUS(FMI2GetEventIndicators(instance, z, nz(S)));
			} else {
				CHECK_STATUS(FMI3GetEventIndicators(instance, z, nz(S)));
			}
		}

		if (isFMI1(S)) {
			nextEventTime = instance->fmi1Functions->eventInfo.nextEventTime;
		} else if (isFMI2(S)) {
			nextEventTime = instance->fmi2Functions->eventInfo.nextEventTime;
		} else {
			nextEventTime = instance->fmi3Functions->nextEventTime;
		}

		z[nz(S)] = nextEventTime - ssGetT(S);
	}
}
#endif


#define MDL_DERIVATIVES
#if defined(MDL_DERIVATIVES)
static void mdlDerivatives(SimStruct *S) {

	logDebug(S, "mdlDerivatives(time=%.16g, majorTimeStep=%d)", ssGetT(S), ssIsMajorTimeStep(S));

	if (isCS(S)) {
		return;  // nothing to do
	}
		
	void **p = ssGetPWork(S);

	FMIInstance *instance = (FMIInstance *)p[0];

    bool inputEvent;

	CHECK_ERROR(setInput(S, true, false, &inputEvent));

    if (isME(S) && inputEvent) {
        ssSetErrorStatus(S, "Unexpected input event in mdlDerivatives().");
    }

	real_T *x = ssGetContStates(S);
	real_T *dx = ssGetdX(S);

	if (isFMI1(S)) {
        CHECK_STATUS(FMI1GetContinuousStates(instance, x, nx(S)));
        CHECK_STATUS(FMI1GetDerivatives(instance, dx, nx(S)));
	} else if (isFMI2(S)) {
        CHECK_STATUS(FMI2GetContinuousStates(instance, x, nx(S)));
        CHECK_STATUS(FMI2GetDerivatives(instance, dx, nx(S)));
	} else {
        CHECK_STATUS(FMI3GetContinuousStates(instance, x, nx(S)));
        CHECK_STATUS(FMI3GetContinuousStateDerivatives(instance, dx, nx(S)));
	}
}
#endif


static void mdlTerminate(SimStruct *S) {

	logDebug(S, "mdlTerminate()");

	void **p = ssGetPWork(S);

	FMIInstance *instance = (FMIInstance *)p[0];

    if (instance) {

	    if (!ssGetErrorStatus(S)) {

		    if (isFMI1(S)) {
		
			    if (isME(S)) {
					if (initialized(S)) {
						CHECK_STATUS(FMI1Terminate(instance));
					}
				    FMI1FreeModelInstance(instance);
			    } else {
					if (initialized(S)) {
						CHECK_STATUS(FMI1TerminateSlave(instance));
					}
				    FMI1FreeSlaveInstance(instance);
			    }

		    } else if (isFMI2(S)) {
				
				if (initialized(S)) {
					CHECK_STATUS(FMI2Terminate(instance));
				}
			    FMI2FreeInstance(instance);
		
		    } else {

				if (initialized(S)) {
					CHECK_STATUS(FMI3Terminate(instance));
				}
			    FMI3FreeInstance(instance);
		
		    }
	    }

	    FMIFreeInstance(instance);
    }

    char *preU = p[3];

    if (preU) {
        free(preU);
    }

	FILE *logFile = (FILE *)p[1];

		if (logFile) {
			fclose(logFile);
			p[1] = NULL;
		}
	}

/*=============================*
* Required S-function trailer *
*=============================*/

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
