/* Do not edit this file. It is produced from the corresponding .m4 source */
/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id: attribute.m4 2873 2017-02-14 02:58:34Z wkliao $ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <pnetcdf.h>
#include <dispatch.h>
#include <pnc_debug.h>
#include <common.h>

/*----< sanity_check_get() >-------------------------------------------------*/
/* This is an independent subroutine. Sanity check for attribute get APIs is
 * simpler, as attribute get APIs are independent subroutines.
 */
static int
sanity_check_get(PNC        *pncp,
                 int         varid,
                 const char *name)
{
    /* check whether variable ID is valid */
    if (varid != NC_GLOBAL && (varid < 0 || varid >= pncp->nvars))
        DEBUG_RETURN_ERROR(NC_ENOTVAR)

    /* sanity check for name */
    if (name == NULL || *name == 0) DEBUG_RETURN_ERROR(NC_EBADNAME)

    if (strlen(name) > NC_MAX_NAME) DEBUG_RETURN_ERROR(NC_EMAXNAME)

    return NC_NOERR;
}

/*----< sanity_check_put() >-------------------------------------------------*/
/* This is a collective subroutine. */
static int
sanity_check_put(PNC        *pncp,
                 int         varid,
                 const char *name,
                 MPI_Offset  nelems,
                 const void *buf)
{
    int err=NC_NOERR;

    /* file should be opened with writable permission */
    if (pncp->flag & NC_MODE_RDONLY)
        DEBUG_RETURN_ERROR(NC_EPERM)

    /* check whether variable ID is valid */
    if (varid != NC_GLOBAL && (varid < 0 || varid >= pncp->nvars))
        DEBUG_RETURN_ERROR(NC_ENOTVAR)

    if (name == NULL || *name == 0) /* name cannot be NULL or NULL string */
        DEBUG_RETURN_ERROR(NC_EBADNAME)

#ifdef NO_NC_GLOBAL_FILLVALUE
    /* See r3403 and RELEASE_NOTES 1.9.0 */
    if (varid == NC_GLOBAL && !strcmp(name, _FillValue))
        DEBUG_RETURN_ERROR(NC_EGLOBAL) /* global _FillValue is not allowed */
#endif

    if (strlen(name) > NC_MAX_NAME) /* name length */
        DEBUG_RETURN_ERROR(NC_EMAXNAME)

    /* check if the name string is legal for netcdf format */
    err = ncmpii_check_name(name, pncp->format);
    if (err != NC_NOERR) return err;

    /* nelems can be zero, i.e. an attribute with only its name */
    if (nelems > 0 && buf == NULL)
        DEBUG_RETURN_ERROR(NC_EINVAL) /* Null arg */

    if (nelems < 0 || (nelems > NC_MAX_INT && pncp->format <= NC_FORMAT_CDF2))
        DEBUG_RETURN_ERROR(NC_EINVAL) /* Invalid nelems */

    return NC_NOERR;
}

/*----< check_EBADTYPE_ECHAR() >---------------------------------------------*/
static int
check_EBADTYPE_ECHAR(PNC *pncp, MPI_Datatype itype, nc_type xtype)
{
    int err;

    /* the max external data type supported by CDF-5 is NC_UINT64 */
    if (xtype <= 0 || xtype > NC_UINT64)
        DEBUG_RETURN_ERROR(NC_EBADTYPE)

    /* For CDF-1 and CDF-2 files, only classic types are allowed. */
    if (pncp->format < NC_FORMAT_CDF5 && xtype > NC_DOUBLE)
        DEBUG_RETURN_ERROR(NC_ESTRICTCDF2)

    /* No character conversions are allowed. */
    err = (((xtype == NC_CHAR) == (itype != MPI_CHAR)) ? NC_ECHAR : NC_NOERR);
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)

    return NC_NOERR;
}

/*----< check_consistency_put() >--------------------------------------------*/
/* This is a collective subroutine and to be called in safe mode. */
static int
check_consistency_put(MPI_Comm      comm,
                      int           varid,
                      const char   *name,
                      nc_type       xtype,
                      MPI_Offset    nelems,
                      const void   *buf,
                      MPI_Datatype  itype,
                      int           err)
{
    int root_name_len, root_varid, minE, rank, mpireturn;
    char *root_name=NULL;
    nc_type root_xtype;
    MPI_Offset root_nelems;

    /* first check the error code, err, across processes */
    TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN, comm);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
    if (minE != NC_NOERR) return minE;

    MPI_Comm_rank(comm, &rank);

    /* check if attribute name is consistent among all processes */
    assert(name != NULL);
    root_name_len = strlen(name) + 1;
    TRACE_COMM(MPI_Bcast)(&root_name_len, 1, MPI_INT, 0, comm);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast root_name_len");

    root_name = (char*) NCI_Malloc((size_t)root_name_len);
    if (rank == 0) strcpy(root_name, name);
    TRACE_COMM(MPI_Bcast)(root_name, root_name_len, MPI_CHAR, 0, comm);
    if (mpireturn != MPI_SUCCESS) {
        NCI_Free(root_name);
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
    }
    if (err == NC_NOERR && strcmp(root_name, name))
        DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_ATTR_NAME)
    NCI_Free(root_name);

    /* check if varid is consistent across all processes */
    root_varid = varid;
    TRACE_COMM(MPI_Bcast)(&root_varid, 1, MPI_INT, 0, comm);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
    if (err == NC_NOERR && root_varid != varid)
        DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_FNC_ARGS)

    /* check if xtype is consistent across all processes */
    root_xtype = xtype;
    TRACE_COMM(MPI_Bcast)(&root_xtype, 1, MPI_INT, 0, comm);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
    if (err == NC_NOERR && root_xtype != xtype)
        DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_ATTR_TYPE)

    /* check if nelems is consistent across all processes */
    root_nelems = nelems;
    TRACE_COMM(MPI_Bcast)(&root_nelems, 1, MPI_OFFSET, 0, comm);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
    if (err == NC_NOERR && root_nelems != nelems)
        DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_ATTR_LEN)

    /* check if buf contents is consistent across all processes */
    if (root_nelems > 0) { /* non-scalar attribute */
        /* note xsz is aligned, thus must use the exact size of buf */
        int itype_size, rank, buf_size;
        void *root_buf;

        MPI_Comm_rank(comm, &rank);
        MPI_Type_size(itype, &itype_size);
        buf_size = (int)root_nelems * itype_size;
        if (rank > 0) root_buf = (void*) NCI_Malloc(buf_size);
        else          root_buf = (void*)buf;

        TRACE_COMM(MPI_Bcast)(root_buf, root_nelems, itype, 0, comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
        if (err == NC_NOERR &&
            (root_nelems != nelems || memcmp(root_buf, buf, buf_size)))
            DEBUG_ASSIGN_ERROR(err, NC_EMULTIDEFINE_ATTR_VAL)
        if (rank > 0) NCI_Free(root_buf);
    }

    /* find min error code across processes */
    TRACE_COMM(MPI_Allreduce)(&err, &minE, 1, MPI_INT, MPI_MIN, comm);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
    if (minE != NC_NOERR) return minE;

    return err;
}





/*----< ncmpi_get_att() >---------------------------------------------------*/
/* This is an independent subroutine.
 * The user buffer data type matches the external type defined in file.
 */
int
ncmpi_get_att(int         ncid,
               int         varid,
               const char *name,
               
                void *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    

    MPI_Datatype itype=MPI_DATATYPE_NULL;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_get(pncp, varid, name);
    if (err != NC_NOERR) return err;

    

    /* calling the subroutine that implements ncmpi_get_att() */
    return pncp->driver->get_att(pncp->ncp, varid, name,
            buf, itype);
}

/*----< ncmpi_get_att_text() >---------------------------------------------------*/
/* This is an independent subroutine.
 * This API never returns NC_ERANGE error, as text is not convertible to numerical types
 */
int
ncmpi_get_att_text(int         ncid,
               int         varid,
               const char *name,
               
                char *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    
    MPI_Datatype itype=MPI_CHAR;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_get(pncp, varid, name);
    if (err != NC_NOERR) return err;

    

    /* calling the subroutine that implements ncmpi_get_att_text() */
    return pncp->driver->get_att(pncp->ncp, varid, name,
            buf, itype);
}

/*----< ncmpi_get_att_schar() >---------------------------------------------------*/
/* This is an independent subroutine.
 *
 */
int
ncmpi_get_att_schar(int         ncid,
               int         varid,
               const char *name,
               
                schar *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    
    MPI_Datatype itype=MPI_SIGNED_CHAR;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_get(pncp, varid, name);
    if (err != NC_NOERR) return err;

    

    /* calling the subroutine that implements ncmpi_get_att_schar() */
    return pncp->driver->get_att(pncp->ncp, varid, name,
            buf, itype);
}

/*----< ncmpi_get_att_uchar() >---------------------------------------------------*/
/* This is an independent subroutine.
 *
 */
int
ncmpi_get_att_uchar(int         ncid,
               int         varid,
               const char *name,
               
                uchar *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    
    MPI_Datatype itype=MPI_UNSIGNED_CHAR;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_get(pncp, varid, name);
    if (err != NC_NOERR) return err;

    

    /* calling the subroutine that implements ncmpi_get_att_uchar() */
    return pncp->driver->get_att(pncp->ncp, varid, name,
            buf, itype);
}

/*----< ncmpi_get_att_short() >---------------------------------------------------*/
/* This is an independent subroutine.
 *
 */
int
ncmpi_get_att_short(int         ncid,
               int         varid,
               const char *name,
               
                short *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    
    MPI_Datatype itype=MPI_SHORT;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_get(pncp, varid, name);
    if (err != NC_NOERR) return err;

    

    /* calling the subroutine that implements ncmpi_get_att_short() */
    return pncp->driver->get_att(pncp->ncp, varid, name,
            buf, itype);
}

/*----< ncmpi_get_att_ushort() >---------------------------------------------------*/
/* This is an independent subroutine.
 *
 */
int
ncmpi_get_att_ushort(int         ncid,
               int         varid,
               const char *name,
               
                ushort *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    
    MPI_Datatype itype=MPI_UNSIGNED_SHORT;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_get(pncp, varid, name);
    if (err != NC_NOERR) return err;

    

    /* calling the subroutine that implements ncmpi_get_att_ushort() */
    return pncp->driver->get_att(pncp->ncp, varid, name,
            buf, itype);
}

/*----< ncmpi_get_att_int() >---------------------------------------------------*/
/* This is an independent subroutine.
 *
 */
int
ncmpi_get_att_int(int         ncid,
               int         varid,
               const char *name,
               
                int *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    
    MPI_Datatype itype=MPI_INT;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_get(pncp, varid, name);
    if (err != NC_NOERR) return err;

    

    /* calling the subroutine that implements ncmpi_get_att_int() */
    return pncp->driver->get_att(pncp->ncp, varid, name,
            buf, itype);
}

/*----< ncmpi_get_att_uint() >---------------------------------------------------*/
/* This is an independent subroutine.
 *
 */
int
ncmpi_get_att_uint(int         ncid,
               int         varid,
               const char *name,
               
                uint *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    
    MPI_Datatype itype=MPI_UNSIGNED;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_get(pncp, varid, name);
    if (err != NC_NOERR) return err;

    

    /* calling the subroutine that implements ncmpi_get_att_uint() */
    return pncp->driver->get_att(pncp->ncp, varid, name,
            buf, itype);
}

/*----< ncmpi_get_att_long() >---------------------------------------------------*/
/* This is an independent subroutine.
 *
 */
int
ncmpi_get_att_long(int         ncid,
               int         varid,
               const char *name,
               
                long *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    
#if SIZEOF_LONG == SIZEOF_INT
    MPI_Datatype itype=MPI_INT;
#elif SIZEOF_LONG == SIZEOF_LONG_LONG
    MPI_Datatype itype=MPI_LONG_LONG_INT;
#endif

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_get(pncp, varid, name);
    if (err != NC_NOERR) return err;

    

    /* calling the subroutine that implements ncmpi_get_att_long() */
    return pncp->driver->get_att(pncp->ncp, varid, name,
            buf, itype);
}

/*----< ncmpi_get_att_float() >---------------------------------------------------*/
/* This is an independent subroutine.
 *
 */
int
ncmpi_get_att_float(int         ncid,
               int         varid,
               const char *name,
               
                float *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    
    MPI_Datatype itype=MPI_FLOAT;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_get(pncp, varid, name);
    if (err != NC_NOERR) return err;

    

    /* calling the subroutine that implements ncmpi_get_att_float() */
    return pncp->driver->get_att(pncp->ncp, varid, name,
            buf, itype);
}

/*----< ncmpi_get_att_double() >---------------------------------------------------*/
/* This is an independent subroutine.
 *
 */
int
ncmpi_get_att_double(int         ncid,
               int         varid,
               const char *name,
               
                double *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    
    MPI_Datatype itype=MPI_DOUBLE;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_get(pncp, varid, name);
    if (err != NC_NOERR) return err;

    

    /* calling the subroutine that implements ncmpi_get_att_double() */
    return pncp->driver->get_att(pncp->ncp, varid, name,
            buf, itype);
}

/*----< ncmpi_get_att_longlong() >---------------------------------------------------*/
/* This is an independent subroutine.
 *
 */
int
ncmpi_get_att_longlong(int         ncid,
               int         varid,
               const char *name,
               
                long long *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    
    MPI_Datatype itype=MPI_LONG_LONG_INT;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_get(pncp, varid, name);
    if (err != NC_NOERR) return err;

    

    /* calling the subroutine that implements ncmpi_get_att_longlong() */
    return pncp->driver->get_att(pncp->ncp, varid, name,
            buf, itype);
}

/*----< ncmpi_get_att_ulonglong() >---------------------------------------------------*/
/* This is an independent subroutine.
 *
 */
int
ncmpi_get_att_ulonglong(int         ncid,
               int         varid,
               const char *name,
               
                unsigned long long *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    
    MPI_Datatype itype=MPI_UNSIGNED_LONG_LONG;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_get(pncp, varid, name);
    if (err != NC_NOERR) return err;

    

    /* calling the subroutine that implements ncmpi_get_att_ulonglong() */
    return pncp->driver->get_att(pncp->ncp, varid, name,
            buf, itype);
}

/*----< ncmpi_put_att() >---------------------------------------------------*/
/* 
 * This is a collective subroutine, all arguments should be consistent among
 * all processes.
 *
 * If attribute name has already existed, it means to overwrite the attribute.
 * In this case, if the new attribute is larger than the old one, this API
 * must be called when the file is in define mode. (This check should be done
 * at the driver.)
 *
 * Note from netCDF user guide:
 * Attributes are always single values or one-dimensional arrays. This works
 * out well for a string, which is a one-dimensional array of ASCII characters.
 *
 * The user buffer data type matches the external type defined in file.
 */
int
ncmpi_put_att(int         ncid,
               int         varid,
               const char *name,
               nc_type xtype,
               MPI_Offset  nelems,   /* number of elements in buf */
               const void *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    

    MPI_Datatype itype=ncmpii_nc2mpitype(xtype);

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_put(pncp, varid, name, nelems, buf);

    /* check NC_EBADTYPE/NC_ECHAR */
    if (err == NC_NOERR) err = check_EBADTYPE_ECHAR(pncp, itype, xtype);

    if (pncp->flag & NC_MODE_SAFE) /* put APIs are collective */
        err = check_consistency_put(pncp->comm, varid, name, xtype, nelems,
                                    buf, itype, err);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_put_att() */
    return pncp->driver->put_att(pncp->ncp, varid, name,
           xtype, nelems, buf, itype);
}

/*----< ncmpi_put_att_text() >---------------------------------------------------*/
/* 
 * This is a collective subroutine, all arguments should be consistent among
 * all processes.
 *
 * If attribute name has already existed, it means to overwrite the attribute.
 * In this case, if the new attribute is larger than the old one, this API
 * must be called when the file is in define mode. (This check should be done
 * at the driver.)
 *
 * Note from netCDF user guide:
 * Attributes are always single values or one-dimensional arrays. This works
 * out well for a string, which is a one-dimensional array of ASCII characters.
 *
 * This API never returns NC_ERANGE error, as text is not convertible to numerical types
 */
int
ncmpi_put_att_text(int         ncid,
               int         varid,
               const char *name,
               
               MPI_Offset  nelems,   /* number of elements in buf */
               const char *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    nc_type xtype=NC_CHAR;
    MPI_Datatype itype=MPI_CHAR;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_put(pncp, varid, name, nelems, buf);

    

    if (pncp->flag & NC_MODE_SAFE) /* put APIs are collective */
        err = check_consistency_put(pncp->comm, varid, name, xtype, nelems,
                                    buf, itype, err);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_put_att_text() */
    return pncp->driver->put_att(pncp->ncp, varid, name,
           xtype, nelems, buf, itype);
}

/*----< ncmpi_put_att_schar() >---------------------------------------------------*/
/* 
 * This is a collective subroutine, all arguments should be consistent among
 * all processes.
 *
 * If attribute name has already existed, it means to overwrite the attribute.
 * In this case, if the new attribute is larger than the old one, this API
 * must be called when the file is in define mode. (This check should be done
 * at the driver.)
 *
 * Note from netCDF user guide:
 * Attributes are always single values or one-dimensional arrays. This works
 * out well for a string, which is a one-dimensional array of ASCII characters.
 *
 *
 */
int
ncmpi_put_att_schar(int         ncid,
               int         varid,
               const char *name,
               nc_type xtype,
               MPI_Offset  nelems,   /* number of elements in buf */
               const schar *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    
    MPI_Datatype itype=MPI_SIGNED_CHAR;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_put(pncp, varid, name, nelems, buf);

    /* check NC_EBADTYPE/NC_ECHAR */
    if (err == NC_NOERR) err = check_EBADTYPE_ECHAR(pncp, itype, xtype);

    if (pncp->flag & NC_MODE_SAFE) /* put APIs are collective */
        err = check_consistency_put(pncp->comm, varid, name, xtype, nelems,
                                    buf, itype, err);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_put_att_schar() */
    return pncp->driver->put_att(pncp->ncp, varid, name,
           xtype, nelems, buf, itype);
}

/*----< ncmpi_put_att_uchar() >---------------------------------------------------*/
/* 
 * This is a collective subroutine, all arguments should be consistent among
 * all processes.
 *
 * If attribute name has already existed, it means to overwrite the attribute.
 * In this case, if the new attribute is larger than the old one, this API
 * must be called when the file is in define mode. (This check should be done
 * at the driver.)
 *
 * Note from netCDF user guide:
 * Attributes are always single values or one-dimensional arrays. This works
 * out well for a string, which is a one-dimensional array of ASCII characters.
 *
 *
 */
int
ncmpi_put_att_uchar(int         ncid,
               int         varid,
               const char *name,
               nc_type xtype,
               MPI_Offset  nelems,   /* number of elements in buf */
               const uchar *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    
    MPI_Datatype itype=MPI_UNSIGNED_CHAR;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_put(pncp, varid, name, nelems, buf);

    /* check NC_EBADTYPE/NC_ECHAR */
    if (err == NC_NOERR) err = check_EBADTYPE_ECHAR(pncp, itype, xtype);

    if (pncp->flag & NC_MODE_SAFE) /* put APIs are collective */
        err = check_consistency_put(pncp->comm, varid, name, xtype, nelems,
                                    buf, itype, err);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_put_att_uchar() */
    return pncp->driver->put_att(pncp->ncp, varid, name,
           xtype, nelems, buf, itype);
}

/*----< ncmpi_put_att_short() >---------------------------------------------------*/
/* 
 * This is a collective subroutine, all arguments should be consistent among
 * all processes.
 *
 * If attribute name has already existed, it means to overwrite the attribute.
 * In this case, if the new attribute is larger than the old one, this API
 * must be called when the file is in define mode. (This check should be done
 * at the driver.)
 *
 * Note from netCDF user guide:
 * Attributes are always single values or one-dimensional arrays. This works
 * out well for a string, which is a one-dimensional array of ASCII characters.
 *
 *
 */
int
ncmpi_put_att_short(int         ncid,
               int         varid,
               const char *name,
               nc_type xtype,
               MPI_Offset  nelems,   /* number of elements in buf */
               const short *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    
    MPI_Datatype itype=MPI_SHORT;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_put(pncp, varid, name, nelems, buf);

    /* check NC_EBADTYPE/NC_ECHAR */
    if (err == NC_NOERR) err = check_EBADTYPE_ECHAR(pncp, itype, xtype);

    if (pncp->flag & NC_MODE_SAFE) /* put APIs are collective */
        err = check_consistency_put(pncp->comm, varid, name, xtype, nelems,
                                    buf, itype, err);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_put_att_short() */
    return pncp->driver->put_att(pncp->ncp, varid, name,
           xtype, nelems, buf, itype);
}

/*----< ncmpi_put_att_ushort() >---------------------------------------------------*/
/* 
 * This is a collective subroutine, all arguments should be consistent among
 * all processes.
 *
 * If attribute name has already existed, it means to overwrite the attribute.
 * In this case, if the new attribute is larger than the old one, this API
 * must be called when the file is in define mode. (This check should be done
 * at the driver.)
 *
 * Note from netCDF user guide:
 * Attributes are always single values or one-dimensional arrays. This works
 * out well for a string, which is a one-dimensional array of ASCII characters.
 *
 *
 */
int
ncmpi_put_att_ushort(int         ncid,
               int         varid,
               const char *name,
               nc_type xtype,
               MPI_Offset  nelems,   /* number of elements in buf */
               const ushort *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    
    MPI_Datatype itype=MPI_UNSIGNED_SHORT;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_put(pncp, varid, name, nelems, buf);

    /* check NC_EBADTYPE/NC_ECHAR */
    if (err == NC_NOERR) err = check_EBADTYPE_ECHAR(pncp, itype, xtype);

    if (pncp->flag & NC_MODE_SAFE) /* put APIs are collective */
        err = check_consistency_put(pncp->comm, varid, name, xtype, nelems,
                                    buf, itype, err);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_put_att_ushort() */
    return pncp->driver->put_att(pncp->ncp, varid, name,
           xtype, nelems, buf, itype);
}

/*----< ncmpi_put_att_int() >---------------------------------------------------*/
/* 
 * This is a collective subroutine, all arguments should be consistent among
 * all processes.
 *
 * If attribute name has already existed, it means to overwrite the attribute.
 * In this case, if the new attribute is larger than the old one, this API
 * must be called when the file is in define mode. (This check should be done
 * at the driver.)
 *
 * Note from netCDF user guide:
 * Attributes are always single values or one-dimensional arrays. This works
 * out well for a string, which is a one-dimensional array of ASCII characters.
 *
 *
 */
int
ncmpi_put_att_int(int         ncid,
               int         varid,
               const char *name,
               nc_type xtype,
               MPI_Offset  nelems,   /* number of elements in buf */
               const int *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    
    MPI_Datatype itype=MPI_INT;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_put(pncp, varid, name, nelems, buf);

    /* check NC_EBADTYPE/NC_ECHAR */
    if (err == NC_NOERR) err = check_EBADTYPE_ECHAR(pncp, itype, xtype);

    if (pncp->flag & NC_MODE_SAFE) /* put APIs are collective */
        err = check_consistency_put(pncp->comm, varid, name, xtype, nelems,
                                    buf, itype, err);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_put_att_int() */
    return pncp->driver->put_att(pncp->ncp, varid, name,
           xtype, nelems, buf, itype);
}

/*----< ncmpi_put_att_uint() >---------------------------------------------------*/
/* 
 * This is a collective subroutine, all arguments should be consistent among
 * all processes.
 *
 * If attribute name has already existed, it means to overwrite the attribute.
 * In this case, if the new attribute is larger than the old one, this API
 * must be called when the file is in define mode. (This check should be done
 * at the driver.)
 *
 * Note from netCDF user guide:
 * Attributes are always single values or one-dimensional arrays. This works
 * out well for a string, which is a one-dimensional array of ASCII characters.
 *
 *
 */
int
ncmpi_put_att_uint(int         ncid,
               int         varid,
               const char *name,
               nc_type xtype,
               MPI_Offset  nelems,   /* number of elements in buf */
               const uint *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    
    MPI_Datatype itype=MPI_UNSIGNED;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_put(pncp, varid, name, nelems, buf);

    /* check NC_EBADTYPE/NC_ECHAR */
    if (err == NC_NOERR) err = check_EBADTYPE_ECHAR(pncp, itype, xtype);

    if (pncp->flag & NC_MODE_SAFE) /* put APIs are collective */
        err = check_consistency_put(pncp->comm, varid, name, xtype, nelems,
                                    buf, itype, err);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_put_att_uint() */
    return pncp->driver->put_att(pncp->ncp, varid, name,
           xtype, nelems, buf, itype);
}

/*----< ncmpi_put_att_long() >---------------------------------------------------*/
/* 
 * This is a collective subroutine, all arguments should be consistent among
 * all processes.
 *
 * If attribute name has already existed, it means to overwrite the attribute.
 * In this case, if the new attribute is larger than the old one, this API
 * must be called when the file is in define mode. (This check should be done
 * at the driver.)
 *
 * Note from netCDF user guide:
 * Attributes are always single values or one-dimensional arrays. This works
 * out well for a string, which is a one-dimensional array of ASCII characters.
 *
 *
 */
int
ncmpi_put_att_long(int         ncid,
               int         varid,
               const char *name,
               nc_type xtype,
               MPI_Offset  nelems,   /* number of elements in buf */
               const long *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    
#if SIZEOF_LONG == SIZEOF_INT
    MPI_Datatype itype=MPI_INT;
#elif SIZEOF_LONG == SIZEOF_LONG_LONG
    MPI_Datatype itype=MPI_LONG_LONG_INT;
#endif

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_put(pncp, varid, name, nelems, buf);

    /* check NC_EBADTYPE/NC_ECHAR */
    if (err == NC_NOERR) err = check_EBADTYPE_ECHAR(pncp, itype, xtype);

    if (pncp->flag & NC_MODE_SAFE) /* put APIs are collective */
        err = check_consistency_put(pncp->comm, varid, name, xtype, nelems,
                                    buf, itype, err);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_put_att_long() */
    return pncp->driver->put_att(pncp->ncp, varid, name,
           xtype, nelems, buf, itype);
}

/*----< ncmpi_put_att_float() >---------------------------------------------------*/
/* 
 * This is a collective subroutine, all arguments should be consistent among
 * all processes.
 *
 * If attribute name has already existed, it means to overwrite the attribute.
 * In this case, if the new attribute is larger than the old one, this API
 * must be called when the file is in define mode. (This check should be done
 * at the driver.)
 *
 * Note from netCDF user guide:
 * Attributes are always single values or one-dimensional arrays. This works
 * out well for a string, which is a one-dimensional array of ASCII characters.
 *
 *
 */
int
ncmpi_put_att_float(int         ncid,
               int         varid,
               const char *name,
               nc_type xtype,
               MPI_Offset  nelems,   /* number of elements in buf */
               const float *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    
    MPI_Datatype itype=MPI_FLOAT;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_put(pncp, varid, name, nelems, buf);

    /* check NC_EBADTYPE/NC_ECHAR */
    if (err == NC_NOERR) err = check_EBADTYPE_ECHAR(pncp, itype, xtype);

    if (pncp->flag & NC_MODE_SAFE) /* put APIs are collective */
        err = check_consistency_put(pncp->comm, varid, name, xtype, nelems,
                                    buf, itype, err);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_put_att_float() */
    return pncp->driver->put_att(pncp->ncp, varid, name,
           xtype, nelems, buf, itype);
}

/*----< ncmpi_put_att_double() >---------------------------------------------------*/
/* 
 * This is a collective subroutine, all arguments should be consistent among
 * all processes.
 *
 * If attribute name has already existed, it means to overwrite the attribute.
 * In this case, if the new attribute is larger than the old one, this API
 * must be called when the file is in define mode. (This check should be done
 * at the driver.)
 *
 * Note from netCDF user guide:
 * Attributes are always single values or one-dimensional arrays. This works
 * out well for a string, which is a one-dimensional array of ASCII characters.
 *
 *
 */
int
ncmpi_put_att_double(int         ncid,
               int         varid,
               const char *name,
               nc_type xtype,
               MPI_Offset  nelems,   /* number of elements in buf */
               const double *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    
    MPI_Datatype itype=MPI_DOUBLE;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_put(pncp, varid, name, nelems, buf);

    /* check NC_EBADTYPE/NC_ECHAR */
    if (err == NC_NOERR) err = check_EBADTYPE_ECHAR(pncp, itype, xtype);

    if (pncp->flag & NC_MODE_SAFE) /* put APIs are collective */
        err = check_consistency_put(pncp->comm, varid, name, xtype, nelems,
                                    buf, itype, err);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_put_att_double() */
    return pncp->driver->put_att(pncp->ncp, varid, name,
           xtype, nelems, buf, itype);
}

/*----< ncmpi_put_att_longlong() >---------------------------------------------------*/
/* 
 * This is a collective subroutine, all arguments should be consistent among
 * all processes.
 *
 * If attribute name has already existed, it means to overwrite the attribute.
 * In this case, if the new attribute is larger than the old one, this API
 * must be called when the file is in define mode. (This check should be done
 * at the driver.)
 *
 * Note from netCDF user guide:
 * Attributes are always single values or one-dimensional arrays. This works
 * out well for a string, which is a one-dimensional array of ASCII characters.
 *
 *
 */
int
ncmpi_put_att_longlong(int         ncid,
               int         varid,
               const char *name,
               nc_type xtype,
               MPI_Offset  nelems,   /* number of elements in buf */
               const long long *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    
    MPI_Datatype itype=MPI_LONG_LONG_INT;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_put(pncp, varid, name, nelems, buf);

    /* check NC_EBADTYPE/NC_ECHAR */
    if (err == NC_NOERR) err = check_EBADTYPE_ECHAR(pncp, itype, xtype);

    if (pncp->flag & NC_MODE_SAFE) /* put APIs are collective */
        err = check_consistency_put(pncp->comm, varid, name, xtype, nelems,
                                    buf, itype, err);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_put_att_longlong() */
    return pncp->driver->put_att(pncp->ncp, varid, name,
           xtype, nelems, buf, itype);
}

/*----< ncmpi_put_att_ulonglong() >---------------------------------------------------*/
/* 
 * This is a collective subroutine, all arguments should be consistent among
 * all processes.
 *
 * If attribute name has already existed, it means to overwrite the attribute.
 * In this case, if the new attribute is larger than the old one, this API
 * must be called when the file is in define mode. (This check should be done
 * at the driver.)
 *
 * Note from netCDF user guide:
 * Attributes are always single values or one-dimensional arrays. This works
 * out well for a string, which is a one-dimensional array of ASCII characters.
 *
 *
 */
int
ncmpi_put_att_ulonglong(int         ncid,
               int         varid,
               const char *name,
               nc_type xtype,
               MPI_Offset  nelems,   /* number of elements in buf */
               const unsigned long long *buf)
{
    int err=NC_NOERR;
    PNC *pncp;
    
    MPI_Datatype itype=MPI_UNSIGNED_LONG_LONG;

    /* check if ncid is valid */
    err = PNC_check_id(ncid, &pncp);
    if (err != NC_NOERR) return err;

    /* sanity check for arguments */
    err = sanity_check_put(pncp, varid, name, nelems, buf);

    /* check NC_EBADTYPE/NC_ECHAR */
    if (err == NC_NOERR) err = check_EBADTYPE_ECHAR(pncp, itype, xtype);

    if (pncp->flag & NC_MODE_SAFE) /* put APIs are collective */
        err = check_consistency_put(pncp->comm, varid, name, xtype, nelems,
                                    buf, itype, err);
    if (err != NC_NOERR) return err;

    /* calling the subroutine that implements ncmpi_put_att_ulonglong() */
    return pncp->driver->put_att(pncp->ncp, varid, name,
           xtype, nelems, buf, itype);
}


