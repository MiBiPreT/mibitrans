import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest
from mibitrans.analysis.to_modflow import MibitransToModflow


class TestMibitransToModflow:
    """Tests conversion from Mibitrans input to FloPy model."""

    @pytest.fixture
    def flopy_object(self, test_hydro_pars, test_att_pars, test_source_pars, test_model_pars):
        """Fixture object of conversion class."""
        mbt_mod = MibitransToModflow(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars)
        mbt_mod.to_modflow()
        return mbt_mod

    @pytest.fixture
    def mf(self, flopy_object):
        """Fixture object of MODFLOW input."""
        return flopy_object.mf

    @pytest.fixture
    def mt(self, flopy_object):
        """Fixture object of MT3DMS input."""
        return flopy_object.mt

    def test_mf_dis(self, mf, test_source_pars, test_model_pars):
        """Test content of MODFLOW dis package."""
        assert mf.dis.nlay == 1
        assert mf.dis.nrow == test_model_pars.model_width * 2 // test_model_pars.dy
        assert mf.dis.ncol == test_model_pars.model_length // test_model_pars.dx
        assert mf.dis.delr[0] == test_model_pars.dx
        assert mf.dis.delc[0] == test_model_pars.dy
        assert mf.dis.top[0, 0] == test_source_pars.depth
        assert mf.dis.botm[0, 0, 0] == 0
        assert mf.dis.perlen[0] == test_model_pars.model_time
        assert mf.dis.nstp[0] == test_model_pars.model_time // test_model_pars.dt

    def test_mf_bas(self, mf, test_hydro_pars, test_source_pars):
        """Test content of MODFLOW bas6 package."""
        assert np.all(mf.bas6.ibound[0, :, 0] == -1)
        assert np.all(mf.bas6.ibound[0, :, -1] == -1)
        assert np.all(mf.bas6.ibound[0, :, 1:-1] == 1)

    def test_mf_lpf(self, mf, test_hydro_pars):
        """Test content of MODFLOW lpf package."""
        assert np.all(mf.lpf.hk[:, :, :] == test_hydro_pars.h_conductivity)

    def test_mt_btn(self, mt, test_hydro_pars):
        """Test content of MODFLOW btn package."""
        mod_src_zone = np.where(mt.btn.sconc[0].array[0, :, 0] > 0)
        assert np.all(mt.btn.icbund[0, mod_src_zone, 0] == -1)
        assert np.all(mt.btn.prsity.array == test_hydro_pars.porosity)

    def test_mt_dsp(self, mt, test_hydro_pars):
        """Test content of MT3DMS dsp package."""
        mbt_al = test_hydro_pars.alpha_x
        assert np.all(mt.dsp.al.array == mbt_al)
        assert np.all(mt.dsp.trpt.array == test_hydro_pars.alpha_y / mbt_al)
        assert np.all(mt.dsp.trpv.array == test_hydro_pars.alpha_z / mbt_al)

    def test_mt_rct(self, mt, test_att_pars):
        """Test content of MT3DMS rct package."""
        assert np.all(mt.rct.rhob.array == test_att_pars.bulk_density)
        assert np.all(
            mt.rct.sp1[0].array == test_att_pars.partition_coefficient * test_att_pars.fraction_organic_carbon
        )
        assert np.all(mt.rct.rc1[0].array == np.float32(test_att_pars.decay_rate))

    def test_run_plot(self, flopy_object):
        """Test if FloPy model runs and if class can produce plot."""
        flopy_object.run_modflow(verbose=False)
        flopy_object.centerline_modflow()
        assert isinstance(plt.gca(), matplotlib.axes._axes.Axes)
