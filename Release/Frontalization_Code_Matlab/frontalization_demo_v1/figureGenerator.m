frontalizers = {'Vito'}
landmarkers = {'dlib' 'ZhuRamanan' 'CMR'}
folder = '/Users/joel/Documents/FaceRepose/FaceYawArray2_front/'
images = dir(strcat(folder, '*.jpg'));
detector = py.dlib.get_frontal_face_detector();
        predictor_file = fullfile(pwd,'shape_predictor_68_face_landmarks.dat')
        predictor = py.dlib.shape_predictor(predictor_file);
        mod = py.importlib.import_module('dlib_detect_script_optimized');
        mod = py.reload(mod)
for i = 1:length(images)
for f = frontalizers
    for l = landmarkers
        name = images(i);
        name = name.name
       makeFrontalFigures(strcat(folder,name), l{1}, f{1})

    end
end
end