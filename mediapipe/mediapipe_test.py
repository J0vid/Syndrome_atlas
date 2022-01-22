
import cv2
import mediapipe as mp
import pandas as pd
import json

mp_face_mesh = mp.solutions.face_mesh

# For static images:
IMAGE_FILES = ['jovid.png']

with mp_face_mesh.FaceMesh(
    static_image_mode=True,
    max_num_faces=1,
    refine_landmarks=True,
    min_detection_confidence=0.5) as face_mesh:
  for idx, file in enumerate(IMAGE_FILES):
    image = cv2.imread(file)
    # Convert the BGR image to RGB before processing.
    results = face_mesh.process(cv2.cvtColor(image, cv2.COLOR_BGR2RGB))


result_lms = results.multi_face_landmarks

keypoints = []
for data_point in results.multi_face_landmarks:
        keypoints.append({
                         'X': data_point(0),
                         'Y': data_point(1),
                         'Z': data_point(2),
                         })

json.dump()


for face_landmarks in results.multi_face_landmarks:
        print('face_landmarks:', face_landmarks)

